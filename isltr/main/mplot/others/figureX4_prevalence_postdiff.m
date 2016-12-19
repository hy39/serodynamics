function [FigH FigL log peak peak_lb peak_ub] = figureX4_prevalence_postdiff()
% 1) Plot the seroconversion fitting the post-wave serology


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
set(FigL, 'Position', [100, 500, 816, 480]);
set(gca,'xtick',[]);
hold;
dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat'); % default titre model
%dat1 = load('out/imm/m4/ph1n1/20160304/mcmc_output_m4.mat'); % default best fit model
%dat1 = load('out/imm/m4_high_inititre/ph1n1/20160315/mcmc_output_m4_final.mat'); % default best fit model
%dat1 = load('out/p0e05/m1.5/ph1n1/20151027/mcmc_output_m1.5_final.mat');
mark = '';
[rmse_1h rmse_1f h1 h2] = plotI(dat1, 1, mark, 1);
[rmse_1h rmse_1f h3] = plotI(dat1, 1 ,mark, 0);

dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat'); % default threshold model
%dat2 = load('out/test/m2/ph1n1/20160313/mcmc_output_m2_final.mat'); % threshold model with cut=1 
%dat2 = load('out/p0e05/m2.2/ph1n1/20151024/mcmc_output_m2.2_final.mat');
%dat2 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat2 = load('out/boost/m3/ph1n1/20160303/mcmc_output_m3_final.mat');
mark = ':';
%[rmse_2h rmse_2f h4] = plotI(dat2, 1, mark, 0);
[rmse_2h rmse_2f h4] = plotI(dat2, 1, mark, 1);
sample_time_K1 = 102;
sample_time_K2 = 235;
line([sample_time_K1 sample_time_K1], [0 30]);
line([sample_time_K2 sample_time_K2], [0 30]);
xlim([92 92+30.5*5]);
ylim([0 30]);
legend([h1 h2 h3 h4],{'Titre model: seroprevalence','Titre model: cumulative incidence','Titre model (no reinfection): cumulative incidence','Threshod model: seroprevalence'})
ylabel('Adjusted seroprevalence & cumulative incidence (%)');

%disp(['rmse (until peak) for mean incidence model1:' num2str(rmse_1h)]);
%disp(['rmse (until peak) for mean incidence model2:' num2str(rmse_2h)]);

%disp(['rmse (full data) for mean incidence model1:' num2str(rmse_1f)]);
%disp(['rmse (full data) for mean incidence model2:' num2str(rmse_2f)]);
function [RMSE_half RMSE_full h1 h2] = plotI(dat, display, marker, mode)
    
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
post = table2array(posteriorTable);
posterior = mean(post);
%posterior = repmat(mean(post),3,1);
%samplesize = 3;
%idx = 1:3;



    %setup initial condition
    if pars.maxi == 2 % only 2 titres
    %  pars.inittitres_flag = 3
    %  pars.inittitres = 0.4;
      [yini age_arr s0_imm] = make_ics_naive2titres( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    else
    %    pars.inittitres_flag = 3
    %    seroconvert_obs = 0.4;
    %    prev = [0.5 0.5 0.5 0.5];
    %    naive = 1 - sum(prev)*seroconvert_obs;
    %    pars.init_prev = [naive prev*seroconvert_obs zeros(1,pars.maxi-length(prev)-1)];
        [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    end
    
    

    
    ab_baseline = Ab.K(init_collect).Abl;
    ab_k = Ab.K(k).Abl;
    for i=1:12
        abid = find(Ab.K(1).numdays>=184+(i-1)*7 & Ab.K(1).numdays<184+i*7)
        ab = ab_baseline(abid);
        ab_age = Ab.K(init_collect).age(abid);
        ser = length(find(ab>=3))./length(ab);
        [yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, ab, ab_age);
        ser1 = gen_seroprev( yini_k1, pars, 1, 1, 4, 0); 
        if ser>0
          ser = ser - 0.0227;
        end
        Seroprev_obs(i) = ser; 
    end
    
    for i=13:60
        %if i == 15
        %    abid = find(Ab.K(2).numdays>=184+(i-1)*14 & Ab.K(2).numdays<184+(i+1)*14)
        %else
        %    abid = find(Ab.K(2).numdays>=184+(i-1)*14 & Ab.K(2).numdays<184+i*14)
        %end
        abid = find(Ab.K(2).numdays>=184+(i-1)*7 & Ab.K(2).numdays<184+i*7)
        ab = Ab.K(2).Abl(abid);
        ab_age = Ab.K(2).age(abid);
        ser = length(find(ab>=3))./length(ab);
        if ser>0
          ser = ser - 0.0227;
        end
        Seroprev_obs(i) = ser; 
        %Seroprev_obs(16) = 0;
    end
    
    
    % recalculate beta to fit R0
    if pars.maxi == 2
    R0 = [];
    vars = posteriorTable.Properties.VariableNames;
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(p));
           %[pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p));
        end
    end
    R0 = calculateR0_fromPars(pars);
    R0_ratio = 1.22./mean(R0);  
    end

    
    vars = posteriorTable.Properties.VariableNames;
   
    for p=1:length(vars)
        if strcmpi('LLH',vars(p))
        else
           [pars] = setParameters(pars,char(vars(p)),posterior(p));
           if pars.maxi == 2
           if strfind(char(vars(p)),'beta') 
               %if posterior(idx(i),p) > 3  
                 %[pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p)*R0_ratio);
                 %R0_new(i) = calculateR0_fromPars(pars);
                  if mode == 0
                      [pars] = setParameters(pars,char(vars(p)),0.0494*1);
                  elseif mode == 1
                      %[pars] = setParameters(pars,char(vars(p)),0.0494*0.976);
                      [pars] = setParameters(pars,char(vars(p)),0.0494*0.979);
                  end
                  %[pars] = setParameters(pars,char(vars(p)),0.0494); % mean from titre model
                  
           end
           if strfind(char(vars(p)),'ContFrac1')
               [pars] = setParameters(pars,char(vars(p)),5.0068);  % mean from titre model
           end
           end
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
    

    
    %[Seroprev] = gen_seroprev( yini_k1, pars, 1, 1, 4, 0); 
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
    javaaddpath('e:\Documents\Github\serodynamics\isltr\java\matlabjava.jar'); %use this one for calculating CI
    import matlabjava.*
    mepar_2 = matlabjava.ParametersSR;
    meser_2 = matlabjava.SerologySR;
    %mepar_2 = matlabjava.Parameters;
    %meser_2 = matlabjava.Serology;
    % set parameters
    meser_2.setParameters(mepar_2);
    meser_2.updateParameters('s0_imm',pars.s0_imm);
    if pars.maxi > 2 & mode == 0
        pars.wan = 0;
    end
    meser_2.updateParameters('wan',pars.wan);
    %meser_2.updateParameters('wan',0);
    meser_2.updateParameters('maxi', pars.maxi);
    meser_2.updateParametersG(pars.arrg);
    meser_2.updateParametersH(pars.arrh);
    meser_2.updateParametersM(pars.matM);
    meser_2.updateParametersBeta(pars.beta);
    
    
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
    %[total_infecteds sero log CI_T1T2 total_age RMSE(i)] = plot_Incidence( yfull(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
    [total_infecteds CInfecteds]= gen_total_infecteds(z, times+1, pars, strain, a);
    [total_infecteds1 CInfecteds1]= gen_total_infecteds(z, times+1, pars, strain, 1);
    [total_infecteds2 CInfecteds2]= gen_total_infecteds(z, times+1, pars, strain, 2);
    [total_infecteds3 CInfecteds3]= gen_total_infecteds(z, times+1, pars, strain, 3);
    [total_infecteds4 CInfecteds4]= gen_total_infecteds(z, times+1, pars, strain, 4);
    if ~exist('peak')
        tmp_peak = find(total_infecteds == max(total_infecteds));
        if length(tmp_peak) == 1
            peak(1) = tmp_peak;
        end
    else
        tmp_peak = find(total_infecteds == max(total_infecteds));
        if length(tmp_peak) == 1
            peak(end+1) = tmp_peak;
        end
    end
   



 if pars.maxi == 2
 % mean(R0_new)
 end
    if samplesize>1
        Y_mean = mean(Y_posterior);
        %Y_mean = reshape(Y_mean(1,:,:),[length(Y_mean(1,:,1)) length(Y_mean(1,1,:))]);
        Y_mean = squeeze(Y_mean);
    else
        Y_mean = Y_posterior;
    end
    Y_lb = quantile(Y_posterior,0.025);
    Y_lb = reshape(Y_lb(1,:,:),[length(Y_lb(1,:,1)) length(Y_lb(1,1,:))]);
    Y_ub = quantile(Y_posterior,0.975);
    Y_ub = reshape(Y_ub(1,:,:),[length(Y_ub(1,:,1)) length(Y_ub(1,1,:))]);

    
%% subplot figb
%figure;
agegroup = 1:4;
T_rel = times(1:365)+1;  %The relative days after 120d 
    
if exist(marker)
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, 0, marker);
    %[total_infecteds sero log CI_T1T2 total_age] = plot_line( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, marker);
    [total_infecteds sero log CI_T1T2 total_age RMSE_half RMSE_full] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);
    
    Serodiff = Seroprev - [0; Seroprev(1:end-1)];
    maxday = find(Serodiff == max(Serodiff));
    tmp_peak_inc = find(total_infecteds == max(total_infecteds));
    tmp_peak_sero = find(Serodiff == max(Serodiff));
    strain = 1;
    %[total_infecteds CInfecteds] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, agegroup); % plot CI
    %[total_infecteds1 CInfecteds1] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 1); % plot CI
    %[total_infecteds2 CInfecteds2] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 2); % plot CI
    %[total_infecteds3 CInfecteds3] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 3); % plot CI
    %[total_infecteds4 CInfecteds4] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 4); % plot CI
    h1 = plot((Seroprev/total_age)*100,'b--');
    plot(maxday,(Seroprev(maxday)/total_age)*100,'bo');
    %plot((CInfecteds/total_age)*100, 'g--');
else
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull]= plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, 0);
    %[total_infecteds sero log CI_T1T2 total_age] = plot_line( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup)
    [total_infecteds sero log CI_T1T2 total_age RMSE_half RMSE_full] = plot_Incidence( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup);

    Serodiff = Seroprev - [0; Seroprev(1:end-1)];
    maxday = find(Serodiff == max(Serodiff));
    tmp_peak_inc = find(total_infecteds == max(total_infecteds));
    tmp_peak_sero = find(Serodiff == max(Serodiff));
    strain = 1;
    [total_infecteds CInfecteds] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, agegroup); % plot CI
    [total_infecteds1 CInfecteds1] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 1); % plot CI
    [total_infecteds2 CInfecteds2] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 2); % plot CI
    [total_infecteds3 CInfecteds3] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 3); % plot CI
    [total_infecteds4 CInfecteds4] = gen_total_infecteds(Y_mean(T0:T0-60+515,:), T_rel, pars, strain, 4); % plot CI
    if mode == 1
        h1 = plot((Seroprev/total_age)*100,'b');
        h2 = plot((CInfecteds/total_age)*100, 'g');
        plot(maxday,(Seroprev(maxday)/total_age)*100,'bo');
        plot(maxday,(CInfecteds(maxday)/total_age)*100,'go');
    end
    if mode == 0
        h1 = plot((CInfecteds/total_age)*100, 'g');
        plot(maxday,(CInfecteds(maxday)/total_age)*100,'go');
    end
    
    %plot(71:7*1:71+7*1*(length(Seroprev_obs)-1),Seroprev_obs,'o');

end

set(gca,'xlim',[1 30.5*10]);
set(gca,'XTick',[1:(365-60-1)/10:365-60]);
line([sample_time_K1 sample_time_K1], [0 100]);
line([sample_time_K2 sample_time_K2], [0 100]);

lastsamplingday = pars.SamplingLastDay + 90;
Xl = [1 lastsamplingday];

months = [
          '   ';
          '   ';
          '   ';
          'Aug';
          'Sep';
          'Oct';
          'Nov';
          'Dec';
          '   ';
          '   ';
          '   '
          ];
     
      
      %% Set Text labels 
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xt = [1:(365-60-1)/10:365-60]+15;
      t = text(Xt,Yl(1)*ones(1,length(Xt)),months(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

      % Remove the default labels
      set(gca,'XTickLabel','')

      % Get the Extent of each text object.  This
      % loop is unavoidable.
      for i = 1:length(t)
        ext(i,:) = get(t(i),'Extent');
      end
      
      % Determine the lowest point.  The X-label will be
      % placed so that the top is aligned with this point.
      LowYPoint = min(ext(:,2));

      % Place the axis label at this point
      XMidPoint = Xl(1)+abs(diff(Xl))/2;
      tl = text(XMidPoint,LowYPoint,'', 'VerticalAlignment','top','HorizontalAlignment','center');

      %set(gca(2),'YTickLabel',[0:0.01:0.02]);
%title(['agegroup' num2str(agegroup)]);
       set(gca, 'FontSize', 11.5);
       
       mTextBox = uicontrol('style','text');
       parentColor = [1 1 1];
       set(mTextBox,'String','Month', 'Position', [20 2   820 20],'FontSize',12,'backgroundcolor',parentColor);
end

end