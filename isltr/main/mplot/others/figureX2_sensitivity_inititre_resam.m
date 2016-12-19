function [FigL] = figureX2_sensitivity_inititre_resam(inittitre)
% 1) Plot the seroconversion between 2 models
% 2) calculate rmse
% plot the mean cumulative incidence
% H0 Standard threshold model
% H1 Titre model
% To do list:
% change this file to do resampling, not using the mean value

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
if exist('FigL')
else
    FigL = figure;
end
%FigL2 = figure;
set(FigL, 'Position', [100, 500, 816, 480]);
%set(gca,'xtick',[]);
hold;
dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');
%dat1 = load('out/imm/m4/ph1n1/20160304/mcmc_output_m4.mat'); % Best fit model
%dat1 = load('out/p0e05/m1.5/ph1n1/20151027/mcmc_output_m1.5_final.mat');
display = 0;

% retrieve a new susceptibility array same as threshold model
AbCut = 3; %3= 1:40
tmp_par = InitParameters(); 
tmp_par.immune_flag = 1;
tmp_par = setParameters(tmp_par,'immune_alpha1',AbCut);
tmp_par = setParameters(tmp_par,'immune_alpha2',AbCut);
tmp_par = setParameters(tmp_par,'immune_alpha3',AbCut);
tmp_par = setParameters(tmp_par,'immune_alpha4',AbCut);
tmp_par.arrh(find(tmp_par.arrh>0.6)) = 1;
tmp_par.arrh(find(tmp_par.arrh<0.6)) = 0;
arrh = tmp_par.arrh;
pars = dat1.par;
pars.arrh_tre = tmp_par.arrh;

plotI(dat1, display, 0, '');
plotI(dat1, display, 1, '');
plotI(dat1, display, 2, '');
%dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');
%dat2 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat2 = load('out/boost/m3/ph1n1/20160303/mcmc_output_m3_final.mat');
%mark = ':';
%plotI(dat2, display, mark, cinf);

function [CInfecteds] = plotI(dat, display, inittitre, marker)
    
posteriorTable = dat.PosteriorSamples;
%pars = dat.par;

%pars.inittitres_flag = 3; %1:Defualt initial immunity; 3:manually defined
%seroprevalence 
if inittitre ~= 0
    inittitres_flag = 3
else
    inittitres_flag = 1;
end
if inittitres_flag == 3
    pars.inittitres_flag = 3;
    
    %otions1 10 times higher
    if inittitre == 2 
        seroconvert_obs = 0.3;
        prev = [0.5 0.5 0.5 0.5];
        naive = 1 - sum(prev)*seroconvert_obs;
        pars.init_prev = [naive prev*seroconvert_obs zeros(1,pars.maxi-length(prev)-1)];
    end
    %otions2 exponential decay
    if inittitre == 1
        geo_ra = 0.282;
        ratio = [geo_ra geo_ra^2 geo_ra^3 geo_ra^4 geo_ra^5 geo_ra^6 geo_ra^7 geo_ra^8 geo_ra^9 geo_ra^10];
        init_prev = ratio./sum(ratio);
        %init_prev(10) = 1-sum(init_prev(1:9));
        pars.init_prev = init_prev;
    end
end
post = table2array(posteriorTable);
%posterior = repmat(mean(post),3,1);
%samplesize = 3;
%idx = 1:3;

%resamples
samplesize = 400;
burnIn = 1000;
total = height(posteriorTable(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);
posterior = post;

samplesize = 1;
idx = 1;
    
%use the posterior mean
    posterior_mean = mean(posterior);
    %b  = 0.0494
    %f1 = 5.0068
    
%calculate R0
    R0 = [];
    for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior_mean(idx(i),p));
           %[pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
        end
    end
    R0(i) = calculateR0_fromPars(pars);
    end
    R0_ratio = 1.22./mean(R0); %%change to 1.22 as target R0  
    %R0_ratio = 1.3101;
if inittitres_flag == 1
    R0_ratio = 1;
end
for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior_mean(idx(i),p));
           %[pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
           if inittitre == 1 || inittitre == 2
             [pars] = setParameters(pars,'beta',posterior_mean(idx(i),1)*R0_ratio);
           end
           %if strfind(char(vars(p)),'immune')
           %    if posterior(idx(i),p) > 3 
           %      [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p)-2);
           %    end
           %end
        end
    end

    %calculate R0
    beta0 = posterior(idx(i),1)*R0_ratio;
    newR0 = calculateR0_fromPars(pars);
    options = optimset('FunValCheck','on')
    x = fzero(@calculateR0_fromBeta, beta0, options, pars);
    if inittitres_flag == 1
    else
        pars.beta = x;
    end
    newR0 = calculateR0_fromPars(pars);
    
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
    %meser_2.updateParametersHtre(pars.arrh_tre); % add h for threshold
    
    
    x0 = yini;  
    x0 = [x0 zeros(1,40) zeros(1,40)];
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
        
    [total_infecteds sero log CI_T1T2 total_age CI_T2 , CI_T1, CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull CITreInfectedsFull] = plot_immune( Y_sample(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display);
    for t=1:6 %titre from 1:10 to 1:320
        
        
        for d = 1:length(CInfecteds)
            if d == 1
                before(1,:) = CIImmInfectedsFull(d,:)./CInfecteds(d);
                after(d,:) = CINaiveInfectedsFull(d,:)./CInfecteds(d);
                threshold(d) = sum(CITreInfectedsFull(d,:))./CInfecteds(d);
                sen_tot(d) = calSensitivity(before(d), after(d), t);   
            else
            %Total
            if d == 360
              disp d;
            end
            before(d,:)=(CIImmInfectedsFull(d,:)-CIImmInfectedsFull(d-1,:))./(CInfecteds(d)-CInfecteds(d-1));    %CIBeforeBoosting
            after(d,:)=(CINaiveInfectedsFull(d,:)-CINaiveInfectedsFull(d-1,:))./(CInfecteds(d)-CInfecteds(d-1)); %CIAfterBoosting
            %After(d)=sum((CINaiveInfectedsFull(d,4:end)-CINaiveInfectedsFull(d-1,4:end)));
            %%%VVVV This is wrong. Should divided by something else than
            %%%CInfecteds
            threshold_full(d,:)=(CITreInfectedsFull(d,:)-CITreInfectedsFull(d-1,:))./(CInfecteds(d)-CInfecteds(d-1)); %CIBeforBoostingThreshold  
            threshold(d)=sum(threshold_full(d,:),2);
            sen_tot(d) = calSensitivity(before(d,:), after(d,:), t);   
            end
        end
        before(isnan(before)) = 0
        %beta_ratio = sen_tot./threshold;
        if t==2 && inittitre == 0
          disp t;
          %figure;
          %calculate sensitivity, precision, and predictability for
          %dichotomy boosting
          cimm = sum(before(:,1:3),2);
          cimm_h = sum(before(:,4:end),2);
          %plot(threshold);
          hold on;
          %plot(cimm);
          
          % sensitivity 
          sen_dich = cimm./(cimm+cimm_h);
          sen_dich(isnan(sen_dich)) = 0;
          %plot(sen_dich);
          sen_dich'*(CInfecteds./sum(CInfecteds))
          
          % precision
          cimm(isnan(cimm)) = 0;
          threshold(isnan(threshold)) = 0;
          prec_dich = cimm./threshold';
          prec_dich(isnan(prec_dich)) = 0;
          prec_dich'*(CInfecteds./sum(CInfecteds))
          %plot(cimm./threshold'); 
          
          % predictability

          %After = sum(after(:,4:end),2);
          Threshold_pos = sum(CITreInfectedsFull,2);
          Titre_all = sum(CINaiveInfectedsFull(:,1:end),2);
          beta_ratio = Threshold_pos./Titre_all;
          beta_ratio(isnan(beta_ratio)) = 0;
          beta_ratio'*(CInfecteds./sum(CInfecteds))
          %plot(beta_ratio);
          %plot(sum(CINaiveInfectedsFull(:,4:end),2),'r')
        end
        Sen_tot(s,t) = (sen_tot*total_infecteds)./sum(total_infecteds);
    end
    end

    if samplesize>1
        Sen_tot_mean = mean(Sen_tot(:,:));
        Sen_lb = quantile(Sen_tot(:,:),0.025);
        Sen_lb = squeeze(Sen_lb);
        Sen_ub = quantile(Sen_tot(:,:),0.975);
        Sen_ub = squeeze(Sen_ub);
    if inittitre == 0
        de=0;
        plot([1:6]-de,Sen_tot_mean,'k.-');
        errorbar([1:6]-de,[Sen_tot_mean],abs(Sen_tot_mean-Sen_lb),abs(Sen_tot_mean-Sen_ub),'k');
    elseif inittitre == 1
        de=0;
        plot([1:6]-de,Sen_tot_mean,'r.-');
        errorbar([1:6]-de,[Sen_tot_mean],abs(Sen_tot_mean-Sen_lb),abs(Sen_tot_mean-Sen_ub),'r');
    elseif inittitre == 2
        de=0;
        plot([1:6]+de,Sen_tot_mean,'b.-');
        errorbar([1:6]+de,[Sen_tot_mean],abs(Sen_tot_mean-Sen_lb),abs(Sen_tot_mean-Sen_ub),'b-');
    end
    else
        Sen_tot_mean = Sen_tot;
        if inittitre == 0
            plot([1:6],Sen_tot_mean,'k.-');
        elseif inittitre == 1
            plot([1:6],Sen_tot_mean,'r.-');
        elseif inittitre == 2
            plot([1:6],Sen_tot_mean,'b.-');
        end
    end
    %errorbar(1:6,[Sen_tot_mean],abs(Sen_tot_mean-Sen_lb),abs(Sen_tot_mean-Sen_ub));
    %%add error bar
    
    
    

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