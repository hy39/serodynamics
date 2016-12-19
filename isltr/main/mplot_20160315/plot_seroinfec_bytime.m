function [FigH FigL log peak peak_lb peak_ub] = plot_seroinfec_bytime()
% 1) Plot the seroconversion given infection by time


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
dat1 = load('out/p0e05/m1/ph1n1/20151024/mcmc_output_m1_final.mat');
%dat1 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat1 = load('out/p0e05/m1.5/ph1n1/20151027/mcmc_output_m1.5_final.mat');
display = 0;
cinf = plotI(dat1, display);

%dat2 = load('out/p0e05/m2/ph1n1/20151024/mcmc_output_m2_final.mat');
dat2 = load('out/p0e05/m2.2/ph1n1/20151024/mcmc_output_m2.2_final.mat');
%dat2 = load('out/imm/m4/ph1n1/20160303/mcmc_output_m4.mat');
%dat2 = load('out/boost/m3/ph1n1/20160303/mcmc_output_m3_final.mat');
mark = ':';
%plotI(dat2, display, mark, cinf);

function [CInfecteds] = plotI(dat, display, marker, cinf)
    
posteriorTable = dat.PosteriorSamples;
pars = dat.par;
pars.inittitres_flag = 3; %more partial protected population
post = table2array(posteriorTable);
%posterior = repmat(mean(post),3,1);
%samplesize = 3;
%idx = 1:3;

%resamples
samplesize = 30;
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
           if strfind(char(vars(p)),'immune')
               if posterior(idx(i),p) > 3 
                 [pars] = setParameters(pars,char(vars(p)),posterior(idx(i),p)-2);
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
    
if exist('marker')
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display, marker);
    %ISP = Seroprev./(CInfecteds);
    %ISP = Seroprev./(cinf);
    
    Seroprev_prev = [0; Seroprev(1:end-1)];
    %cinf_prev = [0; cinf(1:end-1)];
    %ISP = (Seroprev-Seroprev_prev)./(cinf-cinf_prev);
    CInfecteds_prev = [0; CInfecteds(1:end-1)];
    ISP = 1-(Seroprev-Seroprev_prev)./(CInfecteds-CInfecteds_prev);

    %ISP(1:60) = NaN;
    ISP_stable = ISP;
    ISP_stable(1:60) = NaN;
    plot(ISP_stable, 'g');
    plot(ISP(1:60), 'g--');
else
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds Seroprev CINaiveInfecteds CIImmInfecteds CINaiveInfectedsFull CIImmInfectedsFull] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, agegroup, display);
   
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds1 Seroprev1 CINaiveInfecteds1 CIImmInfecteds1 CINaiveInfectedsFull1 CIImmInfectedsFull1] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, 1, display); 
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds2 Seroprev2 CINaiveInfecteds2 CIImmInfecteds2 CINaiveInfectedsFull2 CIImmInfectedsFull2] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, 2, display); 
    
    A=CIImmInfectedsFull(150,:)./CInfecteds(150);
    B=CINaiveInfectedsFull(150,:)./CInfecteds(150);
    bar([A' B'], 1)
    %ISP = Seroprev./CInfecteds;
    %Seroprev = [Seroprev(4:end); 0; 0; 0];
    Seroprev_prev = [0; Seroprev(1:end-1)];
    CInfecteds_prev = [0; CInfecteds(1:end-1)];
    %ISP = (Seroprev-Seroprev_prev)./(CInfecteds-CInfecteds_prev);
    
    ISP = Seroprev./(CInfecteds);
    %ISP(1:60) = NaN;

    ISP_stable = ISP;
    ISP_stable(1:60) = NaN;
    plot(ISP_stable, 'b');
    plot(ISP(1:60), 'b--');
    hold on;
    %plot((CInfecteds-CInfecteds_prev), (Seroprev-Seroprev_prev))
    
    
    ISP = Seroprev./CInfecteds;
    ISP_stable = ISP;
    ISP_stable(1:60) = NaN;
    %plot(ISP_stable, 'g');
    %plot(ISP(1:60), 'g--');

    
    
    
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds1 Seroprev1] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, 1, display);
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds2 Seroprev2] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, 2, display);
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds3 Seroprev3] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, 3, display);
    [total_infecteds sero log CI_T1T2 total_age CI_T2 CI_T1 CInfecteds4 Seroprev4] = plot_immune( Y_mean(T0:T0-60+515,:), pars, T_rel, 1, 4, display);
    %Seroprev1 = [Seroprev1(4:end); 0; 0; 0];
    %Seroprev2 = [Seroprev2(4:end); 0; 0; 0];
    %Seroprev3 = [Seroprev3(4:end); 0; 0; 0];
    %Seroprev4 = [Seroprev4(4:end); 0; 0; 0];
    plot(Seroprev1./CInfecteds1, 'r'); 
    plot(Seroprev2./CInfecteds2, 'g'); 
    plot(Seroprev3./CInfecteds3, 'b'); 
    plot(Seroprev4./CInfecteds4, 'k'); 
end




%new_peak = find(total_infecteds == max(total_infecteds)) + 120;
%formatIn = 'mm/dd/yyyy';
%ddd = datenum({'01/01/2009'},formatIn);
%ddd1 = ddd+new_peak-1;
%ddd2 = ddd+mean(peak)-1;

%peak_lb = quantile(peak,0.025);
%peak_ub = quantile(peak,0.975);
%log = [log '\r\n' 'peak occurs at day '  num2str(mean(peak)) '(' num2str(peak_lb) ', '  num2str(peak_ub) ')'];
%log = [log '\r\n' 'convert to date: ' datestr(ddd2)];
%log = [log '\r\n' 'single_peak from average Incidence '  num2str(mean(new_peak)) '(' datestr(ddd1) ')'];
%disp (log);

set(gca,'xlim',[1 30.5*10]);
set(gca,'XTick',[1:(365-60-1)/10:365-60]);
line([sample_time_K1 sample_time_K1], [0 1]);
line([sample_time_K2 sample_time_K2], [0 1]);

lastsamplingday = pars.SamplingLastDay + 90;
Xl = [1 lastsamplingday];

months = [
          'May';
          'Jun';
          'Jul';
          'Aug';
          'Sep';
          'Oct';
          'Nov';
          'Dec';
          'Jan';
          'Feb';
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