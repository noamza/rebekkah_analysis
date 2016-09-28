function PlotAllFiguresUpdated

%fig1
figure;
n=4;m=4;

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
load('3sets G3MD15PF3 data and results.mat')

load('simulated fano factor means.mat')

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis\simulation examples')
%folder = '\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis\simulation examples';

%Ex=load(sprtntf('%s\simulated rate map example id 170.mat',folder));
%Ex2=load(sprtntf('%s\simulated rate map example2.mat', folder));
Ex=load('simulated rate map example id 170.mat');
Ex2=load('simulated rate map example2.mat');

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')


all_fano_factors=combined_means;

subplot(n,m,[1 2 5 6])
hist(fanos, [0.5:30.5]); hold on;
hist(fanos(fanos<=1),1); box off;
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
xlim([0 31])
xlabel('Fano factors')
ylabel('# of cells') 

t = annotation('textbox');
t.String = ' text';
t.EdgeColor= 'none';

load('simulated fano factor means.mat')
rv=mean(fanos);
cd('C:\\Noam\\Dropbox\\GitTechnion\\rebekkah\\');
PlotShuffleHistogram(n,m,[3 4 7 8], all_fano_factors,rv, 500,'Shuffles','Mean Fano factor', 'C')

t = annotation('textbox');
t.String = ' text';
t.EdgeColor= 'none';

inds=[170 1];
count=0;
for h=1:2 
    
    if h==1
    sm= Ex.simulated_rate_map;
    else
      sm=Ex2.simulated_rate_map;  
    end
    
   subplot(n,m,[9]+(count*4))
   rm=rate_mats_all{inds(h)};
   imagesc(rm);  
   axis off; axis image;
   
   
   text(-0.75,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
   
   
   title(sprintf('%0.1f Hz', mean(rm(:)))); 
   
   subplot(n,m,[10]+(count*4))
   imagesc(gaussian_mats{inds(h)}); 
   axis off; axis image;
   
   subplot(n,m,[11]+(count*4))
   imagesc(sm); 
   axis off; axis image;
   title(sprintf('%0.1f Hz', mean(sm(:)))); 
   
   count=count+1;
end

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 11.6], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 11.6])


%end fig1

%fig 2

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\firing variability corrcoefs')
load('corr coef results half sessions G3MD25PF3.mat')
%best_inds=find(all_corrs>0.7 & all_corrs2>0.7 & twoD_zone_corrs>0.7);

figure;

best_inds=[64 200 213 192];
all_corrs(best_inds(1:2))

for ex_num=1:2
    add=(ex_num-1)*3;
    ImagePlot(ex_num,rate_mats_b,rate_mats_e,peak_rates_b,peak_rates_e,best_inds,add,'A')
end
clearvars

%plot same session results
load('corr coef results of same arenas updated.mat')

best_inds=find(all_corrs>0.6);
best_inds=[best_inds(3) 9];
all_corrs(best_inds)

for ex_num=1:2
    add=(ex_num+2-1)*3;
    ImagePlot(ex_num,all_rm_1,all_rm_2,all_rates_1,all_rates_2,best_inds,add, 'B')
end
clearvars

%plot rescaling session results
load('corr coef results of rescaled arenas.mat')

best_inds=find(all_corrs>0.45);
%best_inds=[best_inds(1) best_inds(3)];
all_corrs(best_inds(1:2))

for ex_num=1:2
    add=(ex_num+4-1)*3;
    ImagePlot(ex_num,all_rm_1,all_rm_2,all_rates_1,all_rates_2,best_inds,add, 'C')
end


set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 30.45], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 30.45])

%end fig 2

% Plot FIGURE 3
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\firing variability corrcoefs')

A=load('shuffled split session firing stability analysis.mat');
B=load('shuffled same areans firing stability analysis.mat');
C=load('shuffled rescaled areans firing stability analysis.mat');

D=load('simutated corr coeff split sesh combined results 1000x.mat');
E=load('simulated same arena stability corr means.mat');
F=load('simulated rescaling arena stability corr means.mat');

Ar=load('half sessions G3MD15PF3 info results.mat','all_corrs2')
Br=load('corr coef results of same arenas COMBINED.mat','all_corrs2');
Cr=load('corr coef results of rescaled arenas.mat','all_corrs2');

PlotThreeShuffleHists(A.mean_corrs2,B.mean_corr2,C.mean_corr2, ...
    nanmean(Ar.all_corrs2), nanmean(Br.all_corrs2),nanmean(Cr.all_corrs2),...
    3010, 'Shuffles', '', '', 'Mean Corr Coeff', 'Shuffles', 'Mean Corr Coeff') 

PlotThreeShuffleHistsTwo(D.combined_means,E.combined_means,F.combined_means, ...
    nanmean(Ar.all_corrs2), nanmean(Br.all_corrs2),nanmean(Cr.all_corrs2),...
    3010, 'Simulations', '', '', 'Mean Corr Coeff', 'Simulations', 'Mean Corr Coeff') 

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 17.4])
% end Fig 3

%fig 9:cluster score shuffles
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
load('cluster score shuffled results.mat','mean_cluster_scores_pt_35')
hyper=load('cluster scores.mat', 'cluster_scores_35');
minimum=load('cluster scores minimum field.mat', 'cluster_scores_35');
second=load('cluster scores second max.mat', 'cluster_scores_35');

figure;
n=2;m=2;

subplot(n,m,1);
hist(hyper.cluster_scores_35,[0.1,0.3,0.5,0.7,0.9,1.1])
xlabel('Cluster scores')
ylabel('# of cells') 
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

PlotShuffleHistogram(n,m,2, mean_cluster_scores_pt_35,mean(hyper.cluster_scores_35), 3200,'Shuffles','Mean cluster score', 'B')

a=hyper.cluster_scores_35;
b=second.cluster_scores_35;
c=minimum.cluster_scores_35;
means=[mean(a), mean(b),mean(c)];
stderrs=[std(a)/sqrt(length(a)),std(a)/sqrt(length(a)),std(b)/sqrt(length(c))];

subplot(n,m,3)
errorbar(means,stderrs,'o')
box off;

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 11.6], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 11.6])
%end fig 9

%S5


function ImagePlot(ex_num,rate_mats_b,rate_mats_e,peak_rates_b,peak_rates_e,best_inds,add,letter)

n=6;
m=3;

ind=best_inds(ex_num);

subplot(n,m,1+add)
imagesc(rate_mats_b{ind});
axis off;
axis image;
hold on;

c={{[1,0,0]},{[1,1,0]},{[0,1,0]},{[0,0,1]},{[0,1,1]},{[1,0,1]}};
for h=1:len(mi)
plot(mi(h,2),mi(h,1),'o','MarkerSize',10,'LineWidth',2,'MarkerFaceColor',c);
end

if add+1==1 || add+1==7 || add+1==13
text(-0.45,1.1,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
end

subplot(n,m,2+add)
imagesc(rate_mats_e{ind})
axis off;
   axis image

if ex_num==1
    subplot(n,m,1+add)
    title('First half','HorizontalAlignment', 'right')
    
    subplot(n,m,2+add)
    title('Second half','HorizontalAlignment', 'right')
end

rates_1=peak_rates_b{ind};
[rates_1,rank]=sort(rates_1);
rates_2=peak_rates_e{ind};
rates_2=rates_2(rank);

r=corrcoef(rates_1,rates_2);

max_lim=max([rates_1(end) rates_2(end)]);

subplot(n,m,3+add)
scatter(rates_1,rates_2,[],'k');hold on;
%xlabel('1st half')
ylabel('firing rate [Hz]')
lsline
%ylim([0 max_lim])
%xlim([1 length(rates_1)])
box off;


disp('')




