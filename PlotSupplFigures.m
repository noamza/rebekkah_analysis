function PlotSupplFigures

% S5: border score thresholds

figure;
n=1;
m=2;

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('shuffled max peak distributions.mat', ...
    'max_less_pt_05', 'max_less_pt_15') 

load('3sets G3MD25 data and results.mat', 'norm_dist')

rv=sum(norm_dist<0.05)/length(norm_dist);
PlotShuffleHistogram(n,m,1, max_less_pt_05,rv, 3500,'Shuffles','% cells near border (<0.05)', 'A')
rv=sum(norm_dist<0.15)/length(norm_dist);
PlotShuffleHistogram(n,m,2, max_less_pt_15,rv, 3100,'','% cells near border (<0.15)', 'B')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 8.52], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 8.52])


%S7
figure;
n=1;
m=2;

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\cluster score results')
load('cluster score shuffled results.mat', ...
    'mean_cluster_scores_pt_25','mean_cluster_scores_pt_45') 

load('cluster scores.mat', 'cluster_scores_25','cluster_scores_45')

rv=mean(cluster_scores_25);
PlotShuffleHistogram(n,m,1, mean_cluster_scores_pt_25,rv, 3500,'Shuffles','Mean Cluster score (near=0.25)', 'A')
rv=mean(cluster_scores_45);
PlotShuffleHistogram(n,m,2, mean_cluster_scores_pt_45,rv, 3100,'','Mean Cluster score(near=0.45)', 'B')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 8.52], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 8.52])