function figure6a( PosteriorSamples, pars )
%UNTITLED Summary of this function goes here
%Plot seroconversion date using paired sera
setISL;
global proj Antibody;
init_collect = 1;
second_collect = 2;
third_collect = 3;
k = 2;
Ab = Antibody;

pairs = [TitresTablePaired.Age_1 TitresTablePaired.Levels_T1 TitresTablePaired.Levels_T2];
    
for i=1:pars.maxa
  lage = pars.ages(i,1);
  uage = pars.ages(i,2);
  idx = find(pairs(:,1)<uage & pairs(:,1)>lage);
  pairs_a = pairs(idx,:);
  seroconv_a = find((pairs_a(:,3)-pairs_a(:,2))>=2);
  seroconv_imm_a = find((pairs_a(:,3)-pairs_a(:,2))>=2 & pairs_a(:,2)>=1);
  no_total_a(i) = length(pairs_a(:,1));
  ratio_sero_a(i) = length(seroconv_a) / no_total_a(i);
  ratio_sero_imm_a(i) = length(seroconv_imm_a) / no_total_a(i); 
end
  seroconv = find((pairs(:,3)-pairs(:,2))>=2);
  seroconv_imm = find((pairs(:,3)-pairs(:,2))>=2 & pairs(:,2)>=1);
  no_total = length(pairs(:,1));
  ratio_sero = length(seroconv) / no_total;
  ratio_sero_imm = length(seroconv_imm) / no_total;
  
  x_mean = [ratio_sero ratio_sero_a];
  y_mean = [ratio_sero_imm ratio_sero_imm_a];
  
    H = figure;
    set(H, 'Position', [500, 500, 820, 540]);
    %x - y
    H = bar([y_mean' x_mean'-y_mean'],'stacked', 'BarWidth',0.5,'LineWidth',0.1);
    myC = [255/255 0 0; 60/255 160/255 230/255];
    for k=1:2
        set(H(k),'facecolor',myC(k,:))
    end
    hold on;    
    %errorbar([1:5],y_mean,abs(xerrlb), xerrub,'k.','Marker', 'none');
    %errorbar([1:5],x_mean,abs(yerrlb), yerrub,'k.','Marker', 'none');
    

    set(gca,'XTickLabel',{'Total';'3-19';'20-39';'40-64';'65+'});
    ylabel('Proportion');
    xlabel('Age groups');
end

