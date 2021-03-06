function [FigH FigL log peak peak_lb peak_ub] = figure3_sensitivity_age()
% 1) Plot the seroconversion between 2 models
% 2) calculate rmse
% plot the mean cumulative incidence
% H0 Standard threshold model
% H1 Titre model


p = path;
%path(p,'../');
path(p,'lib/');

global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;

%retrieve parameters from posterior
%model:3
if exist('samplesize') == 0
    samplesize = 10;
end
if exist('burnIn') == 0
    burnIn = 1000;
end

%MLE approach
%idx = find(posterior(:,llhidx)==max(posterior(:,llhidx)));
%samplesize = length(idx);






%Mean parameters
FigL = figure;
%FigL2 = figure;
set(FigL, 'Position', [100, 500, 816, 480]);
%set(gca,'xtick',[]);
hold;
%dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');
dat1 = load('out/imm/m4/ph1n1/20160304/mcmc_output_m4.mat'); % Best fit model
%dat1 = load('out/p0e05/m1.5/ph1n1/20151027/mcmc_output_m1.5_final.mat');
display = 0;
cinf = plotI(dat1, display);

%dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');
%dat2 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat2 = load('out/boost/m3/ph1n1/20160303/mcmc_output_m3_final.mat');
%mark = ':';
%plotI(dat2, display, mark, cinf);

function [CInfecteds] = plotI(dat, display, marker, cinf)
    
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
%pars.inittitres_flag = 3; %1:Defualt initial immunity; 3:more partial protected population
post = table2array(posteriorTable);
%posterior = repmat(mean(post),3,1);
%samplesize = 3;
%idx = 1:3;

%resamples
samplesize = 10;
burnIn = 1000;
total = height(posteriorTable(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);
posterior = post;

for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
           %if strfind(char(vars(p)),'immune')
           %    if posterior(idx(i),p) > 3 
           %      [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p)-2);
           %    end
           %end
        end
    end

    %set parameters
    beta = pars.beta;
    AbB = [pars.AbB1 pars.AbB2 pars.AbB3 pars.AbB4];
    immune_alpha = [pars.immune_alpha1 pars.immune_alpha2 pars.immune_alpha3 pars.immune_alpha4];
    lastsamplingday = pars.SamplingLastDay + 90;% -60 -> +30
    if sum(pars.arrh(1,1:4,2))<0.0000001
      pars.arrh(1,1:4,2);
    end
    
    %setup initial condition
    ab_baseline = Ab.K(init_collect).Abl;
    ab_k = Ab.K(k).Abl;
    if pars.maxi == 2 % only 2 titres
        [yini age_arr s0_imm] = make_ics_naive2titres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
        ab_baseline(find(ab_baseline <3)) = 0;
        ab_baseline(find(ab_baseline >2)) = 1;
        ab_k(find(ab_k <3)) = 0;
        ab_k(find(ab_k >2)) = 1;
    else
        [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    end

    %setep simulation time
    T0 = pars.OutbreakStartingDay;
    meanKdays(1) = mean(pars.Antibody.K(1).numdays - T0);
    meanKdays(2) = mean(pars.Antibody.K(2).numdays - T0);
    sample_time_K1 = round(meanKdays(1));
    sample_time_K2 = round(meanKdays(2));
    times = 0:1:lastsamplingday;
    sample_size_K1 = Ab.K(1).samplesize;
    sample_size_K2 = Ab.K(2).samplesize;

    %run simulation
    %initialize objects
    %javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar; %use jave working dir 
    %javaaddpath(pars.javapath); %set ./java as default dir
    javaaddpath('e:\Documents\Github\serodynamics\isltr\java\matlabjava.jar');
    import matlabjava.*
    mepar_2 = matlabjava.ParametersSR;
    meser_2 = matlabjava.SerologySR;
    % set parameters
    meser_2.setParameters(mepar_2);
    meser_2.updateParameters('s0_imm',pars.s0_imm);
    meser_2.updateParameters('wan',pars.wan);
    %meser_2.updateParameters('wan',0);
    meser_2.updateParameters('maxi', pars.maxi);
    meser_2.updateParametersG(pars.arrg);
    meser_2.updateParametersH(pars.arrh);
    meser_2.updateParametersM(pars.matM);
    meser_2.updateParametersBeta(pars.beta);
    
    
    x0 = yini;  
    %x0 = [x0 zeros(1,40)];
    x0 = [x0 zeros(1,80)];
    %create an array starting from 1 Jan
    yfull = zeros(T0+lastsamplingday,length(x0(1,:))); %make a full y array storing simulated data from day 1 in the year
    [t y] = ode23(@(t,x)odef_islmodjava(t,x, meser_2), times, x0);  
    clear('mepar_2');
    clear('meser_2');
    NDA = 0;
    if isfield(pars,'OutbreakNDA') == 1
      NDA = pars.OutbreakNDA;
    end
    
    yfull(1:T0,:) = repmat(x0,T0,1); %create initial data until T0
    if NDA < 0 
        yfull(1:T0-NDA,:) = repmat(x0,T0-NDA,1); %create initial data until T0-NDA
    end
    
    yfull(T0-NDA:T0-NDA+length(y(:,1))-1,:) = y; %save output into the array 
    T = t;
    a=1:pars.maxa;
    Y_posterior(i,:,:) = yfull; 
    
    %calculate peak date
    z = squeeze(yfull);
    strain = 1;
    agegroup = 1:4;
    T_rel = times(1:365)+1; 

   

end
    %% Figure1 Sensitivity of total age group for the full model with 3 different initial titres
    %%
    for s=1:samplesize
        Y_sample = Y_posterior(s,:,:);
        Y_sample = squeeze(Y_sample);
        
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull] = plot_immune( Y_sample(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display);
    for t=1:6 %titre from 1:10 to 1:320
        for d = 10:length(CInfecteds)
            %Total
            before=(CIImmInfectedsFull(d,:)-CIImmInfectedsFull(d-1,:))./(CInfecteds(d)-CInfecteds(d-1));    %CIBefore
            after=(CINaiveInfectedsFull(d,:)-CINaiveInfectedsFull(d-1,:))./(CInfecteds(d)-CInfecteds(d-1)); 
            sen_tot(d) = calSensitivity(before, after, t);   
        end
        Sen_tot(s,t) = (sen_tot*total_infecteds)./sum(total_infecteds);
    end
        Sen_tot_mean = mean(Sen_tot(:,:));
        Sen_lb = quantile(Sen_tot(:,:),0.025);
        Sen_lb = squeeze(Sen_lb);
        Sen_ub = quantile(Sen_tot(:,:),0.975);
        Sen_ub = squeeze(Sen_ub);
    end
    plot(Sen_tot_mean,'o-');
    %%add error bar
    
    
    %% Figure2 Sensitivity of each age group for the full model with the default initial titres
    for s=1:samplesize
    for a=1:4
        
        [total_infecteds_a sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds_a Seroprev_a CINaiveInfecteds_a CIImmInfecteds_a CINaiveInfectedsFull_a CIImmInfectedsFull_a] = plot_immune( Y_sample(T0:T0-60+515,:), pars, T_rel, 1, a, display); 
        for t=1:6 %titre from 1:10 to 1:320
            for d = 10:length(CInfecteds)
                 %Age
                 before_a=(CIImmInfectedsFull_a(d,:)-CIImmInfectedsFull_a(d-1,:))./(CInfecteds_a(d)-CInfecteds_a(d-1));    %CIBefore
                 after_a=(CINaiveInfectedsFull_a(d,:)-CINaiveInfectedsFull_a(d-1,:))./(CInfecteds_a(d)-CInfecteds_a(d-1)); 
                 sen(a,d) = calSensitivity(before_a, after_a, t);
            end
            Sen(s,a,t) = (sen(a,:)*total_infecteds_a)./sum(total_infecteds_a);
        end
    end
    end
    
    
    for a=1:4
        if samplesize>1
            Sen_mean = mean(Sen(:,a,:));
            %Y_mean = reshape(Y_mean(1,:,:),[length(Y_mean(1,:,1)) length(Y_mean(1,1,:))]);
            Sen_mean = squeeze(Sen_mean);
        else
            Sen_mean = Sen(a,:);
        end
        Sen_bar(:,a) = Sen_mean;
    end
    b = bar(Sen_bar);
    b(1).FaceColor = 'r';
    b(2).FaceColor = 'g';
    b(3).FaceColor = 'b';
    b(4).FaceColor = 'y';
    
    color_code = {'k','k','k','k'};
    for a=1:4
        if samplesize>1
            Sen_mean = mean(Sen(:,a,:));
            Sen_mean = squeeze(Sen_mean);
        else
            Sen_mean = Sen(a,:);
        end
        Sen_lb = quantile(Sen(:,a,:),0.025);
        Sen_lb = squeeze(Sen_lb);
        Sen_ub = quantile(Sen(:,a,:),0.975);
        Sen_ub = squeeze(Sen_ub);
        
        errorbar([1+(a-2.5)*0.18:1:6+(a-2.5)*0.18],[Sen_mean],abs(Sen_mean-Sen_lb),abs(Sen_mean-Sen_ub),color_code{a});
    end

    disp ('done');
end






function [sensitivity] = calSensitivity(CIBefore, CIAfter, AbCut)
    cutid = AbCut + 1; 
    X_i = CIBefore;
    Y_i = CIAfter;
    X = sum(X_i(cutid:end)); %seropositivity
    Y = sum(Y_i(cutid:end)); %seropositivity
    TP = (1-X)*((Y-X)./(1-X));
    %FN = X;
    sensitivity = TP;
      
end





end