function [] = calDIC()
% Summary of figur2
% Plot disease and serological dynamics
% For each age group, 5 subfigures are plotted in a row.
% fig2a, heat map of immune dynamics
% fig2b, disease dynamics
% no resampling but onlty mle
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot all age group information.
% Return [FigH FigL log peak peak_lb peak_ub]
% also return the likelihood ratio test
% H0 Standard threshold model
% H1 Titre model
% Actual cumulative incidence is defined as the average of CI(titre)&CI(threshold)

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

%Mean parameters
%dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');
%dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');
dat1 = load('out/p0e05/m1.12/ph1n1/20151024/mcmc_output_m1.12_final.mat');
dat2 = load('out/p0e05/m2.12/ph1n1/20151024/mcmc_output_m2.12_final.mat');


function [] = plotI(dat, display, marker)
    
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
post = table2array(posteriorTable);
%resamples
samplesize = 20;
burnIn = 1000;
total = height(posteriorTable(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);
posterior = post;


    Antibody = par.Antibody;
    SampleSize = Antibody.samplesize;
    Abl = Antibody.Abl;
    age = Antibody.age;
    numdays = Antibody.numdays;
    Abl(find(Abl>par.maxt)) = par.maxt; %substitute Ab level >maxt

for i = 1:samplesize
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
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
    [yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, ab_baseline, Ab.K(init_collect).age);
    [yini_k2 age_arr_k2] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, ab_k, Ab.K(k).age);

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
    javaaddpath(pars.javapath); %set ./java as default dir
    import matlabjava.*
    mepar_2 = matlabjava.Parameters;
    meser_2 = matlabjava.Serology;
    % set parameters
    meser_2.setParameters(mepar_2);
    meser_2.updateParameters('s0_imm',pars.s0_imm);
    meser_2.updateParameters('wan',pars.wan);
    meser_2.updateParameters('maxi', pars.maxi);
    meser_2.updateParametersG(pars.arrg);
    meser_2.updateParametersH(pars.arrh);
    meser_2.updateParametersM(pars.matM);
    meser_2.updateParametersBeta(pars.beta);
    
    
    x0 = yini;  
    %create an array starting from 1 Jan
    yfull = zeros(T0+lastsamplingday,length(x0(1,:))); %make a full y array storing simulated data from day 1 in the year
    [t xt] = ode23(@(t,x)odef_islmodjava(t,x, meser_2), times, x0);  
    clear('mepar_2');
    clear('meser_2');

    T_rel = (0:1:SamplingLastDay)+1; %%The relative days after 120d 
    for a = 1:par.maxa
        serodist = gen_strain_titres(Xt, T_rel, a, par); % [times x titre]
        serodist(find(serodist< 0)) = 0; %replace negative values by zero
        multinomial(a).p = serodist./(diag(sum(serodist,2))*ones(size(serodist)));
        multi_p_byagebytime(a,:,:) = multinomial(a).p;
    end
    
        for a=1:par.maxa %for each age group
        prob_x = multinomial(a).p;
        prob = prob_x * par.matY;
        obs_titres = yt_array(a).obs_titres;
        llh_array = log(prob).*obs_titres;
        llh_array(isinf(llh_array)) = -10E2; %If -Inf, replace to -10E2 
        llh_array(~isfinite(llh_array))=0;  %If NaN, replace to zero
        LLH_Map_age(a).llh = sum(sum(llh_array));
        end
    LLH_Map_totalage.llh = sum([LLH_Map_age.llh]);
    nLLH = -LLH_Map_totalage.llh;
    
       
 
end

end
end