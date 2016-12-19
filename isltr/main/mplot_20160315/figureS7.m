function [] = figureS7()
%Summary of this function goes here
%Adapted from calR0_byOutput 
%theta_beta: average beta (redundant, think whether I need it or not)
%par: from model output
%mean_posteior: average value of posteriorsamples
%This function won't produce credible interval
%1 Feb 2015
% Written by Sean Yuan (hyuan@imperial.ac.uk)

%chech here for filling the area between the curves,
%http://stackoverflow.com/questions/6245626/matlab-filling-in-the-area-between-two-sets-of-data-lines-in-one-figure
%http://stackoverflow.com/questions/19910982/matlab-fill-area-between-lines

figure;

%%the first model
out_dir1 = ['out/m4.101/ph1n1/20150721/'];
load([out_dir1 'Rt_m4.101.mat']);
Rt_lower = quantile(Rt,0.025);
Rt_upper = quantile(Rt,0.975);
Rt_mean = mean(Rt)
%subplot(3,1,2);
plot(Rt_mean,'k');
hold on;
baseLine = 0.2;
x = 1:360;
h1 = fill([x fliplr(x)], [Rt_upper fliplr(Rt_lower)], 'b','EdgeColor','none');
alpha(.5);

out_dir4 = ['out/m4.41/ph1n1/20150806/'];
load([out_dir4 'Rt_m4.41.mat']);
Rt_lower = quantile(Rt,0.025);
Rt_upper = quantile(Rt,0.975);
Rt_mean = mean(Rt)
plot(Rt_mean,'k');
baseLine = 0.2;
x = 1:360;
h2 = fill([x fliplr(x)], [Rt_upper fliplr(Rt_lower)], 'r','EdgeColor','none');
alpha(.5);
addLabel();
%% add the ticks 
legend([h1 h2],{'\hat{R}_{B}','R_{B}(15% protection)'});  % Only the blue and green lines appear


end

function [] = addLabel() 
%% add the ticks 
legend({'$\bar{R}_{B}$','$R_{B}$ CI','$\bar{R}_{SIR}$','$R_{SIR}$ CI','$\bar{R}_{SIR}$','$R_{SIR}$ CI'},'Interpreter','latex')

lastsamplingday = 305;
xlim([1 lastsamplingday]);
ylim([0.6 1.4]);
set(gca,'XTick',[1:30.5:lastsamplingday]);
%set('YTickLabel',);

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
     
      ax = axis;    % Current axis limits
      axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
      Yl = ax(3:4);  % Y-axis limits
      % Place the text labels
      Xt = [1:(365-60-1)/10:365-60]+15;
      t = text(Xt,Yl(1)*ones(1,length(Xt)),months(1:1:11,:));
      set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
      
      % Remove the default labels
      set(gca,'XTickLabel','')
end
