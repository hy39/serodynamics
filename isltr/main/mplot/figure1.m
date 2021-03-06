function [FigH chi_P1 chi_P2] = figure1( PosteriorSamples, pars, burnIn, samplesize)
% Plot sera histogram at T1 and T2 with uncertainty for each age groups
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
% Plot each age group information for only 10 titres levels
% E:\working\Projects.IC\Projects\isl\mat\Misltr\isltr-1.4\out\mcmc\ph1n1\20150106

global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;
posterior = table2array(PosteriorSamples);
%retrieve parameters from posterior
%model:3
if exist('samplesize') == 0
    samplesize = 10;
end
if exist('burnIn') == 0
    burnIn = 1000;
end
post = mean(posterior(burnIn:end,:));
total = length(posterior(:,1))-burnIn;
idx = burnIn + round(rand(1, samplesize) * total);

for i = 1:length(idx)
    vars = PosteriorSamples.Properties.VariableNames;
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
    lastsamplingday = pars.SamplingLastDay;

    %setup initial condition
    [yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
    [yini_k1 age_arr_k1] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(init_collect).Abl, Ab.K(init_collect).age);
    [yini_k2 age_arr_k2] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(k).Abl, Ab.K(k).age);
    [yini_k3 age_arr_k3] = make_ics_fromtitres_byage( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu, Ab.K(k).Abl, Ab.K(k).age);


    %setep simulation time
    T0 = pars.OutbreakStartingDay;
    meanKdays(1) = mean(pars.Antibody.K(1).numdays - T0);
    meanKdays(2) = mean(pars.Antibody.K(2).numdays - T0);
    meanKdays(3) = mean(pars.Antibody.K(3).numdays - T0);
    sample_time_K1 = round(meanKdays(1));
    sample_time_K2 = round(meanKdays(2));
    times = 0:1:lastsamplingday;
    sample_size_K1 = Ab.K(1).samplesize;
    sample_size_K2 = Ab.K(2).samplesize;

    %run simulation
    %initialize objects

    %javaaddpath e:\workspace\MyJavaProject\bin\matlabjava.jar; %use jave working dir 
    %javaaddpath(pars.javapath); %set ./java as default dir
    javaaddpath java\matlabjava.jar;
    mepar_3b = matlabjava.Parameters;
    meser_3b = matlabjava.Serology;
    meser_3b.setParameters(mepar_3b);
    % set parameters
    meser_3b.updateParametersG(pars.arrg);
    meser_3b.updateParametersH(pars.arrh);
    meser_3b.updateParametersM(pars.matM);
    meser_3b.updateParametersBeta(pars.beta);  
    meser_3b.updateParameters('wan',pars.wan);
    meser_3b.updateParameters('s0_imm', pars.s0_imm);
    x0 = yini;  
    [t y] = ode23(@(t,x)odef_islmodjava(t,x, meser_3b), times, x0);  
    %[t y] = ode23(@(t,x)odef_islmod(t,x,pars), times, x0);
    clear('mepar_3b');
    clear('meser_3b');
    T = t;
    %a=1:pars.maxa;
    for a=1:4
    Xout_k1(:,:) = retrieve_histogram(y, pars, times, sample_time_K1, a); % model output
    Xoutput_k1_list(a,i,:) = Xout_k1;
    Xout_k2 = retrieve_histogram(y, pars, times, sample_time_K2, a); % model output
    Xoutput_k2_list(a,i,:) = Xout_k2;
    end
    
    

    atot=1:pars.maxa;
    Xout_k1_tot = retrieve_histogram(y, pars, times, sample_time_K1, atot); % model output
    Xoutput_k1_list_tot(i,:) = Xout_k1_tot;
    Xout_k2_tot = retrieve_histogram(y, pars, times, sample_time_K2, atot); % model output
    Xoutput_k2_list_tot(i,:) = Xout_k2_tot;

end

FigH = figure;
set(FigH, 'Position', [150, 150, 980, 1220]);
hold on;
plot_all_sera();
plot_age_sera();

%% plot age specific sera 
function [] = plot_age_sera()
%plot disease dynamics
%FigH = figure;
%set(FigH, 'Position', [100, 100, 1200, 960]);


%http://uk.mathworks.com/matlabcentral/newsreader/view_thread/158636
age_label = {'<20','20-39','40-64','\geq65'};





for a=1:4;
strain = 1;
%%%VVV Still working 20150109
% Calculate error rates
xlb = reshape(quantile(Xoutput_k1_list(a,:,:),0.05),1,pars.maxi); %column vector repersents the posterior samples
xub = reshape(quantile(Xoutput_k1_list(a,:,:),0.95),1,pars.maxi);
xlb_k2 = reshape(quantile(Xoutput_k2_list(a,:,:),0.05),1,pars.maxi); %column vector repersents the posterior samples
xub_k2 = reshape(quantile(Xoutput_k2_list(a,:,:),0.95),1,pars.maxi);


%% Plot 1st column
% susceptible population dynamics from model output
subplot(5,2,(a)*2+1);



%Xout_k1 = retrieve_histogram(y, pars, times, sample_time_K1, a); % model output
Xout_k1 = reshape(mean(Xoutput_k1_list(a,:,:)),1,pars.maxi);
Yout_k1 = retrieve_histogram(yini_k1, pars, 0, 1, a); % observed data at day 0
Yout3 = [Yout_k1 Xout_k1']; %Observed / Model output

%chi_P1 = table2_goodnessfit(Yout_k1', Xout_k1',sample_size_K1, sample_size_K1); 

titres = 0:1:pars.maxi; % x axis
count_titres = Yout3(1:pars.maxi,:); % y axis
nf = 5;
count_titres(1,:) = count_titres(1,:)/nf; % normalize titres level 0
for i=pars.maxi:-1:2
  count_titres(i+1,:) = count_titres(i,:); 
end
count_titres(2,:) = 0;

[hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');

set(hax(1),'YLim',[0 1/nf + 0.000001]);
set(hax(1),'YTick',[0:1/(nf*5):1/nf]);
set(hax(1),'YTickLabel',[0:20:100]);

set(hax(2),'YLim',[0 1/nf + 0.000001]);
set(hax(2),'YTick',[0:1/(nf*5):1/nf]);
set(hax(2),'YTickLabel',[0:100/(nf*5):100/nf]);
set(hax(2),'YTickLabel',[]);

set(hax(1),'XLim',[-0.5 10.5]);
set(hax(2),'XLim',[-0.5 10.5]);
set(hax(1),'XTick',[0:10]);
set(hax(1),'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ']);
set(hax(2), 'Ticklength', [0 0]);
set(hbar1,'FaceColor', [0.5 0.5 0.5], 'EdgeColor', ['none'], 'barWidth', 0.3);
set(hbar2,'FaceColor', [0.1 0.4 0.8], 'EdgeColor', ['none'], 'barWidth', 0.3);
%if a == 1
%    title(['All']);
%else
    title([char(age_label(a)) '']);
%end
if a==2
    axes(hax(1)); ylabel('Percentage of undetectable titre','FontSize', 14);
else
    axes(hax(1)); ylabel('');
end
axes(hax(2)); ylabel(''); % Disable the ylabel
xlabel('Titres', 'FontSize', 12);

%legend('Observed','Model output');

%%%VVV Still working 20150109
%%%Don't plot the line and marker at mean
xerrlb = xlb-Xout_k1;
xerrub = xub-Xout_k1;
hold(hax(1), 'on');
errorbar(hax(1), titres+0.3-0.15, count_titres(:,2), -[xerrlb(1)/5 0 xerrlb(2:end)], [xerrub(1)/5 0 xerrub(2:end)],'k.','Marker', 'none');
%hold(hax(2), 'on');
%errorbar(hax(2), titres-0.15, count_titres(:,1), -[err(1) 0 err(2:end)], -[err(1) 0 err(2:end)],'g.','Marker', 'none');

% 
%error bar for Yout_k1 using binomial distribution
options = optimoptions('fsolve','Display','off'); % Option to display output
for i=1:length(Yout_k1)
    pbin = Yout_k1(i)/sum(Yout_k1);
    y0 = pbin;
    sizeage = length(find(Ab.K(1).age>=pars.ages(a,1) & Ab.K(1).age<pars.ages(a,2)));
    if pbin>0
        x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    else
        x = 0;
    end
    pbinlb(i) = x;
    yerrlb(i) = pbinlb(i) - pbin;
    % find pmax
    x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    pbinub(i) = x;
    yerrub(i) = pbinub(i) - pbin;
end

hold(hax(2), 'on');
errorbar(hax(2), titres-0.15, count_titres(:,1), -[yerrlb(1)/5 0 yerrlb(2:end)], [yerrub(1)/5 0 yerrub(2:end)],'k.','Marker', 'none');
Titres = [
          '<1:10 ';
          '      ';
          '1:10  ';
          '1:20  ';
          '1:40  ';
          '1:80  ';
          '1:160 ';
          '1:320 ';
          '1:640 ';
          '1:1280';
          '1:2560'
          ];
      
        %% Set Text labels 
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xl = [1:11];
      Xt = [0.3:10.3];
      t = text(Xt,Yl(1)*ones(1,length(Xt)),Titres(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'FontSize',9);
      set(gca,'XTickLabel','')

%% Plot 2nd column
% susceptible population dynamics from model output
subplot(5,2,(a)*2+2);
%Xout_k2 = retrieve_histogram(y, pars, times, sample_time_K2, a); % model output
%ylim([0 0.1]);
title('model output');
Xout_k2 = reshape(mean(Xoutput_k2_list(a,:,:)),1,pars.maxi);
Yout_k2 = retrieve_histogram(yini_k2, pars, 0, 1, a);
Yout4 = [Yout_k2 Xout_k2'];

%chi_P2 = table2_goodnessfit(Yout_k2', Xout_k2',sample_size_K2, sample_size_K2); 
titres = 0:1:pars.maxi; % x axis
count_titres = Yout4(1:pars.maxi,:); % y axis
ratio = 0.1;
count_titres(1,:) = count_titres(1,:)/nf; % normalize titres level 0
for i=pars.maxi:-1:2
  count_titres(i+1,:) = count_titres(i,:); 
end
count_titres(2,:) = 0;

[hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');
set(hax(1),'YLim',[0 1/nf + 0.000001]);
set(hax(1),'YTick',[0:1/(nf*5):1/nf]);
set(hax(1),'YTickLabel',[0:20:100]);
set(hax(1),'YTickLabel',[]);

set(hax(2),'YLim',[0 1/nf + 0.000001]);
set(hax(2),'YTick',[0:1/(nf*5):1/nf]);
set(hax(2),'YTickLabel',[0:100/(nf*5):100/nf]);

set(hax(1),'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ']);
set(hax(1),'XLim',[-0.5 10.5]);
set(hax(2),'XLim',[-0.5 10.5]);
set(hax(1),'XTick',[0:10]);
set(hax(2), 'Ticklength', [0 0]);
set(hbar1,'FaceColor', [0.5 0.5 0.5], 'EdgeColor', ['none'], 'barWidth', 0.3);
set(hbar2,'FaceColor', [0.1 0.4 0.8], 'EdgeColor', ['none'], 'barWidth', 0.3);
axes(hax(1)); ylabel(''); %Disable the labels
title([char(age_label(a)) ' yrs']);
if a==2
    axes(hax(2)); ylabel('Percentage of detectable titres', 'FontSize', 14);
else
    axes(hax(2)); ylabel('');
end
xlabel('Titres');

xerrlb = xlb_k2-Xout_k2;
xerrub = xub_k2-Xout_k2;
hold(hax(1), 'on');
errorbar(hax(1), titres+0.3-0.15, count_titres(:,2), -[xerrlb(1)/5 0 xerrlb(2:end)], [xerrub(1)/5 0 xerrub(2:end)],'k.','Marker', 'none');


options = optimoptions('fsolve','Display','off'); % Option to display output
for i=1:length(Yout_k2)
    pbin = Yout_k2(i)/sum(Yout_k2);
    y0 = pbin;
    % find pmin
    sizeage = length(find(Ab.K(2).age>=pars.ages(a,1) & Ab.K(2).age<pars.ages(a,2)));
    if pbin>0
        x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    else
        x = 0;
    end
    pbinlb(i) = x;
    yerrlb(i) = pbinlb(i) - pbin;
    % find pmax
    x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    pbinub(i) = x;
    yerrub(i) = pbinub(i) - pbin;
end

hold(hax(2), 'on');
errorbar(hax(2), titres-0.15, count_titres(:,1), -[yerrlb(1)/5 0 yerrlb(2:end)], [yerrub(1)/5 0 yerrub(2:end)],'k.','Marker', 'none');
Titres = [
          '<1:10 ';
          '      ';
          '1:10  ';
          '1:20  ';
          '1:40  ';
          '1:80  ';
          '1:160 ';
          '1:320 ';
          '1:640 ';
          '1:1280';
          '1:2560'
          ];
      
      %% Set Text labels 
      ax = axis;     % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xl = [1:11];
      Xt = [0.3:10.3];
      t = text(Xt,Yl(1)*ones(1,length(Xt)),Titres(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'FontSize',9);
      set(gca,'XTickLabel','')
if a==4
   %text(0.5, 0.02, 'test', ...
   %'VerticalAlignment', 'bottom', ...
   %'HorizontalAlignment', 'center');
   
mTextBox = uicontrol('style','text');
   parentColor = [1 1 1];
   %set(mTextBox,'String',[age_label{4}], 'Position', [210 60   160 30],'FontSize',12,'backgroundcolor',parentColor);
   %set(mTextBox,'String','Titres (Baseline)', 'Position', [210 60   160 30],'FontSize',12,'backgroundcolor',parentColor);
   mTextBox2 = uicontrol('style','text');
   set(mTextBox2,'String','Titres', 'Position', [420 60   160 30],'FontSize',12,'backgroundcolor',parentColor);
end

end

end

%% plot total sera 
function [] = plot_all_sera
    
xlb = quantile(Xoutput_k1_list_tot,0.025); %column vector repersents the posterior samples
xub = quantile(Xoutput_k1_list_tot,0.975);
xlb_k2 = quantile(Xoutput_k2_list_tot,0.025); %column vector repersents the posterior samples
xub_k2 = quantile(Xoutput_k2_list_tot,0.975);


a=1:pars.maxa;
strain = 1;

% Plot 1st column
% susceptible population dynamics from model output
subplot(5,2,1);


Xout_k1 = mean(Xoutput_k1_list_tot);
Yout_k1 = retrieve_histogram(yini_k1, pars, 0, 1, a); % observed data at day 0
Yout3 = [Yout_k1' Xout_k1']; %Observed / Model output

chi_P1 = table2_goodnessfit(Yout_k1', Xout_k1',sample_size_K1, sample_size_K1); 

titres = 0:1:pars.maxi; % x axis
count_titres = Yout3(1:pars.maxi,:); % y axis
nf = 5;
count_titres(1,:) = count_titres(1,:)/nf; % normalize titres level 0
for i=pars.maxi:-1:2
  count_titres(i+1,:) = count_titres(i,:); 
end
count_titres(2,:) = 0;

[hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');


set(hax(1),'YLim',[0 1/nf + 0.000001]);
set(hax(1),'YTick',[0:1/(nf*5):1/nf]);
set(hax(1),'YTickLabel',[0:20:100]);

set(hax(2),'YLim',[0 1/nf + 0.000001]);
set(hax(2),'YTick',[0:1/(nf*5):1/nf]);
set(hax(2),'YTickLabel',[0:100/(nf*5):100/nf]);
set(hax(2),'YTickLabel',[]);
set(hax(1),'XLim',[-0.5 10.5]);
set(hax(2),'XLim',[-0.5 10.5]);
set(hax(1),'XTick',[0:10]);
set(hax(1),'XTickLabel',[' ';' '; ' ';' ';' ';' ';' ';' ';' ';' ';' ']);
set(hax(2), 'Ticklength', [0 0]);
set(hbar1,'FaceColor', [0.5 0.5 0.5], 'EdgeColor', ['none'], 'barWidth', 0.3);
set(hbar2,'FaceColor', [0.1 0.4 0.8], 'EdgeColor', ['none'], 'barWidth', 0.3); %[0.1 0.1 0.5] -> [0.1 0.4 0.8]
title(['\fontsize{14}Baseline']);
axes(hax(1)); ylabel('');
axes(hax(2)); ylabel('');
%xlabel('Titres', 'FontSize', 12);
legend('Observed','Model output');

%%%VVV Still working 20150109
%%%Don't plot the line and marker at mean
xerrlb = xlb-Xout_k1;
xerrub = xub-Xout_k1;
hold(hax(1), 'on');
errorbar(hax(1), titres+0.3-0.15, count_titres(:,2), -[xerrlb(1)/10 0 xerrlb(2:end)], [xerrub(1)/10 0 xerrub(2:end)],'k.','Marker', 'none');



%error bar for Yout_k1 using binomial distribution
options = optimoptions('fsolve','Display','off'); % Option to display output
for i=1:length(Yout_k1)
    pbin = Yout_k1(i)/sum(Yout_k1);
    y0 = pbin;
    % find pmin
    sizeage = sample_size_K1;
    if pbin>0
        x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    else
        x = 0;
    end
    pbinlb(i) = x;
    yerrlb(i) = pbinlb(i) - pbin;
    % find pmax
    x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    pbinub(i) = x;
    yerrub(i) = pbinub(i) - pbin;

end

hold(hax(2), 'on');
errorbar(hax(2), titres-0.15, count_titres(:,1), -[yerrlb(1)/10 0 yerrlb(2:end)], [yerrub(1)/10 0 yerrub(2:end)],'k.','Marker', 'none');

Obs_titres.k1 = [count_titres(:,1) count_titres(:,1)+[yerrlb(1)/10 0 yerrlb(2:end)]' count_titres(:,1)+[yerrub(1)/10 0 yerrub(2:end)]'];

Titres = [
          '<1:10 ';
          '      ';
          '1:10  ';
          '1:20  ';
          '1:40  ';
          '1:80  ';
          '1:160 ';
          '1:320 ';
          '1:640 ';
          '1:1280';
          '1:2560'
          ];
      
        %% Set Text labels 
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xl = [1 11];
      Xt = [0.3:10.3];
      t = text(Xt,Yl(1)*ones(1,length(Xt)),Titres(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'FontSize',9);
      % Remove the default labels
      set(gca,'XTickLabel','')


% Plot 2nd column
% susceptible population dynamics from model output
subplot(5,2,2);
%Xout_k2 = retrieve_histogram(y, pars, times, sample_time_K2, a); % model output
%ylim([0 0.1]);
title('model output');
Xout_k2 = mean(Xoutput_k2_list_tot);
Yout_k2 = retrieve_histogram(yini_k2, pars, 0, 1, a);
Yout4 = [Yout_k2' Xout_k2'];

chi_P2 = table2_goodnessfit(Yout_k2', Xout_k2',sample_size_K2, sample_size_K2); 

titres = 0:1:pars.maxi; % x axis
count_titres = Yout4(1:pars.maxi,:); % y axis
ratio = 0.1;
count_titres(1,:) = count_titres(1,:)/nf; % nStepslize titres level 0
for i=pars.maxi:-1:2
  count_titres(i+1,:) = count_titres(i,:); 
end
count_titres(2,:) = 0;

[hax,hbar1,hbar2] = plotyy(titres+0.3-0.15,count_titres(:,2),titres-0.15,count_titres(:,1),'bar','bar');

set(hax(1),'YLim',[0 1/nf + 0.000001]);
set(hax(1),'YTick',[0:1/(nf*5):1/nf]);
set(hax(1),'YTickLabel',[0:20:100]);
set(hax(1),'YTickLabel',[]);

set(hax(2),'YLim',[0 1/nf + 0.000001]);
set(hax(2),'YTick',[0:1/(nf*5):1/nf]);
set(hax(2),'YTickLabel',[0:100/(nf*5):100/nf]);


set(hax(1),'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ']);
set(hax(1),'XLim',[-0.5 10.5]);
set(hax(2),'XLim',[-0.5 10.5]);
set(hax(1),'XTick',[0:10]);
set(hax(2), 'Ticklength', [0 0]);
set(hbar1,'FaceColor', [0.5 0.5 0.5], 'EdgeColor', ['none'], 'barWidth', 0.3);
set(hbar2,'FaceColor', [0.1 0.4 0.8], 'EdgeColor', ['none'], 'barWidth', 0.3);
%title(['Baseline']);
axes(hax(1)); ylabel('');
axes(hax(2)); ylabel('');
title('\fontsize{14}Follow-up');
legend('Observed','Model output');

xerrlb = xlb_k2-Xout_k2;
xerrub = xub_k2-Xout_k2;
hold(hax(1), 'on');
errorbar(hax(1), titres+0.3-0.15, count_titres(:,2), -[xerrlb(1)/10 0 xerrlb(2:end)], [xerrub(1)/10 0 xerrub(2:end)],'k.','Marker', 'none');



for i=1:length(Yout_k2)
    pbin = Yout_k2(i)/sum(Yout_k2);
    y0 = pbin;
    sizeage = sample_size_K2;
    x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    pbinlb(i) = x;
    yerrlb(i) = pbinlb(i) - pbin;
    % find pmax
    x = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    pbinub(i) = x;
    yerrub(i) = pbinub(i) - pbin;
end

hold(hax(2), 'on');
errorbar(hax(2), titres-0.15, count_titres(:,1), -[yerrlb(1)/10 0 yerrlb(2:end)], [yerrub(1)/10 0 yerrub(2:end)],'k.','Marker', 'none');
Obs_titres.k2 = [count_titres(:,1) count_titres(:,1)+[yerrlb(1)/10 0 yerrlb(2:end)]' count_titres(:,1)+[yerrub(1)/10 0 yerrub(2:end)]'];

Titres = [
          '<1:10 ';
          '      ';
          '1:10  ';
          '1:20  ';
          '1:40  ';
          '1:80  ';
          '1:160 ';
          '1:320 ';
          '1:640 ';
          '1:1280';
          '1:2560'
          ];
      
        %% Set Text labels 
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xl = [1:11];
      Xt = [0.3:10.3];
      t = text(Xt,Yl(1)*ones(1,length(Xt)),Titres(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'FontSize',9);
      set(gca,'XTickLabel','')    
    
end
        
        

end