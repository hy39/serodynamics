function [ ] = plot_boosting_chart()
%Plot antibody boosting probability
ab_rate = [2 3 4];
max_titre = 9;
for i=1:length(ab_rate)
    pa.AbB = ab_rate(i);
    Boosting(i,:) = [poisspdf([0:max_titre-1],pa.AbB)/sum(poisspdf([0:max_titre-1],pa.AbB))];
end
bar(Boosting',1);
set(gca,'XTickLabel',{'1:10','1:20','1:40','1:80','1:160','1:320','1:640','1:1280','1:2560'});


x = 1:max_titre;
y = 1:max_titre; 
X = linspace(min(x),max(x)+1,length(x)+1); % return the index of x on 2D grid
Y = linspace(min(y),max(y)+1,length(y)+1); % return the index of y on 2D grid
%y = [0.1 0.5 3 3.2 3.5 6 7 8 9];
[X,Y] = meshgrid([x max(x)+1],[y max(y)+1]);
boosting = zeros(max(x),max(y));

for j=1:max_titre
boosting(j:max(y),j) = flipud(Boosting(j).abl);
end
%add a row
boosting(max(y)+1,:) = boosting(max(y),:);
boosting(:,max(x)+1) = boosting(:,max(x));

h = surf(X,Y,boosting);
view(0,90);
set(gca,'xlim',[1 11]);
set(gca,'XTickLabel',[0:10]);
set(gca,'ylim',[1 11]);
set(gca,'YTickLabel',[0:10]);
end