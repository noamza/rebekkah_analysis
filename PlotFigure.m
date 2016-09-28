function PlotFigure

% Plot FIGURE 3
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\firing variability corrcoefs')

A=load('shuffled split session firing stability analysis.mat');
B=load('shuffled same session firing stability analysis.mat');
C=load('shuffled rescaled session firing stability analysis.mat');

Ar=load('corr coef results half sessions G3MD25PF3.mat');
Br=load('corr coef results of same arenas updated.mat');
Cr=load('corr coef results of rescaled arenas.mat');

PlotThreeShuffleHists(A.mean_corrs2,B.mean_corr2,C.mean_corr2, ...
    nanmean(Ar.all_corrs2), nanmean(Br.all_corrs2),nanmean(Cr.all_corrs2),...
    3300, 'Shuffles', '', '', 'Mean Corr Coeff', 'Shuffles', 'Mean Corr Coeff') 
% end Fig 3

%Figure 6
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat','norm_dist')
load('shuffled max peak distributions.mat', 'max_less_pt_1')

figure;
n=3;m=2;

subplot(n,m,2)
PlotHistogram(norm_dist, 'Normalized distance','Number of cells','B')

subplot(n,m,3)
PlotShuffleHistogram(max_less_pt_1, '% of cells','Shuffles', 'C',sum(norm_dist<0.1)/length(norm_dist),3000)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 17.4])

clearvars

% Figure 7
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat','norm_second_dist')
load('shuffled max peak distributions.mat', 'second_max_less_pt_1')
load('border vs center mean rates ns wo hyper.mat')
%load('border conjunctivity shuffled.mat','means_30')

figure;
n=2;m=4;

subplot(n,m,[1 2])
PlotShuffleHistogram(second_max_less_pt_1, 'Normalized distance','Number of cells',...
    'A',sum(norm_second_dist<0.1)/length(norm_second_dist),3000)

subplot(n,m,3)
x=border_mean;
stderrb=nanstd(x)/sqrt(length(~isnan(x)));
x=central_mean;
stderrc=nanstd(x)/sqrt(length(~isnan(x)));
errorbar([nanmean(border_mean) nanmean(central_mean)], [stderrb stderrc],'o') 
hold on;
errorbar([nanmean(border_mean)], stderrb,'o') 
ylabel('Mean firing rate')
box off;
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

subplot(n,m,[5 6])
scatter(border_mean,central_mean,'o'); hold on;
xlabel('Mean rate [Hz] (border)')
ylabel('Mean rate [Hz] (central)')
axis equal; box off;
x=0:80;
xlim([0 80])
ylim([0 80])
plot(x,x,'-')
text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)


% subplot(n,m,[6 7])
%  PlotShuffleHistogram(border_over_total_30, '% of cells','Shuffles',...
%      'D',mean(means_30),3000)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 11.6], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 11.6])



function PlotHistogram(data, x_label, y_label,letter)
hist(data);
box off;
xlabel(x_label)
ylabel(y_label)
text(-0.25,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

function PlotShuffleHistogram(data, x_label, y_label,letter,real_value,y_lim)
hist(data); hold on;
y=0:y_lim;
x=ones(1,y_lim+1)*real_value;
plot(x,y,'-')
xlabel(x_label)
ylabel(y_label)
ylim([0 y_lim])
box off;
text(-0.25,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)