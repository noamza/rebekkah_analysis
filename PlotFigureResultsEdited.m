% figure 1: fanofactors

cd('C:\Users\Dori\Desktop\sorted rates') %opens examples 
parms.dir_load_data = 'C:\Users\Dori\Desktop\sorted rates';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

figure;
n= 2;
m= 4;

figure_num= [{'A'} {'B'} {'C'} {'D'} {'E'} {'F'} {'G'} {'H'}]; 

for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
   
    order= 1:length(S.sorted_means);
    
    subplot(n,m,i);
    plot(order, S.sorted_means, 'ko-');
    title(sprintf('Fano factor= %0.2f', std(S.sorted_means)/mean(S.sorted_means)));
    box off;
    set(gca,'ygrid','on')
    
    if i==1 || i==5
       ylabel('firing rate'); 
    end
  
    text(-0.35,1.05, figure_num(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    
end

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 11.6], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 11.6])

%end figure 1: fano factors


%load('border vs center mean rates.mat')
%load('border fields in top max firing smoothed data.mat')
%load('border conjunctivity shuffled smoothed data.mat')

%figure 2
cd('N:\users\rebekkah\results and info of analysis')

load('variability and border distances.mat')
load('simulations.mat')
load('shuffled max peak distributions.mat')

figure;
n=3;
m=2;

subplot(n,m,1); %fano factors histogram
hist(fano_factor, 20); hold on; box off;
xlabel('Fano factors');
ylabel('Number of cells');
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
xlim([0 27])
ylim([0 33])

subplot(n,m,2); %hyperfield distances from border histogram
hist(norm_dist); hold on; box off;
xlabel('Normalized distance');
text(-0.25,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
xlim([0 0.5])
ylim([0 50])

subplot(n,m,3); %fano factor simulation data
hist(fanos_med); hold on; box off;
y= 0:150;
x= ones(1,151)* median(fano_factor); 
plot(x,y,'LineWidth', 2, 'Color', [0.85 0.16 0]);
ylabel('Simulations');
xlabel('Median Fano factor');
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
xlim([0.1 2.4])
ylim([0 140])

subplot(n,m,4); %distances shuffled data
hist(max_less_pt_1); hold on; box off;
y= 0:300;
x= ones(1,301)* sum(norm_dist<0.1)/length(norm_dist); 
plot(x,y,'-', 'LineWidth', 2, 'Color', [0.85 0.16 0]);
ylabel('Shuffles')
xlabel('Percentage of cells')
text(-0.25,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
xlim([0 0.7])
ylim([0 270])

subplot(n,m,6); %distances simulated data
hist(sums_ns); hold on; box off;
y= 0:200;
x= ones(1,201)* sum(norm_dist<0.1)/length(norm_dist); 
plot(x,y,'-','LineWidth', 2, 'Color', [0.85 0.16 0]);
xlabel('Percentage of cells')
ylabel('Simulations')
text(-0.25,1.05,'F','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
xlim([0.2 0.7])
ylim([0 130])

cd('N:\users\rebekkah\results and info of analysis')
load('norm variability analysis results.mat')

subplot(n,m,5)

hold on;
box off;
xlabel('Normalized field rank')
ylabel('Normalized firing rate')
 

all_norm_orders= all_norm_orders';
all_sorted_rates= all_sorted_rates';
f= fit(all_norm_orders,all_sorted_rates,'poly1');
p1= 0.7748;
p2= 0.0674;
p1_1= 0.7423;
p2_1= 0.04754;
p1_2= 0.8074;
p2_2= 0.08726;

x=0:0.1:1.1;

y=(p1*x +p2);
y1= (p1_1*x +p2_1);
y2= (p1_2*x +p2_2);

h= plot(x,y, 'Color', [0.85 0.16 0], 'LineWidth', 2); hold on;
h1= plot(x,y1,'--', 'Color', [0.85 0.16 0],'LineWidth', 2);
h2= plot(x,y2,'--','Color', [0.85 0.16 0],'LineWidth', 2);
text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

d=errorbar([0.5:0.1:1], means_by_bins, stderr, 'k.');

xlim([0 1.1])
ylim([0 1.1])
%legend([d h h1],{'mean firing rates','linear fit', '95% conf. intervals'});

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 17.4])
% end figure 2

% Figure S3 relating to 2C (shuffled max peak)
load('shuffled max peak distributions.mat')
load('variability and border distances.mat')

figure;
n=1;
m=2;

subplot(n,m,1)
hist(max_less_pt_05); hold on; box off;
y= 0:280;
x= ones(1,281)* sum(norm_dist < 0.05)/length(norm_dist);
plot(x,y,'-', 'LineWidth', 2);
xlabel('percentage of cells near border (<0.05)');
ylabel('shuffles');
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

subplot(n,m,2)
hist(max_less_pt_15, 8); hold on; box off;
y=0:290;
x= ones(1,291)* sum(norm_dist < 0.15)/length(norm_dist);
plot(x,y,'-', 'LineWidth', 2);
xlabel('percentage of cells near border (<0.15)');
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)




set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 5.8], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 5.8])

%end of Figure S3

%Figure 3
load('border fields in top max firing hyperfield excluded.mat')
load('border conjunctivity shuffled UPDATED hyperfield excluded.mat')
load('border vs center mean rates UPDATED.mat')
load('shuffled max peak distributions.mat')

figure;
n=2;
m=2;

subplot(n,m,1)
hist(second_max_less_pt_1); hold on; box off;
y= 0:300;
x= ones(1,301)* 0.5581; 
plot(x,y,'-', 'LineWidth', 2);
ylabel('shuffles')
xlabel('percentage of cells')
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

subplot(n,m,2)
plot(border_mean, central_mean, 'ko'); hold on;
axis equal; box off;
x=0:40;
plot(x,x,'k-', 'LineWidth', 2);
xlabel('mean firing rate (border)');
ylabel('mean firing rate (center)');
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

subplot(n,m,3)
hist(means_25); hold on; box off;
y= 0:300;
x= ones(1,301)* mean(border_over_total_25);
plot(x,y,'-', 'LineWidth', 2);
ylabel('shuffles');
xlabel('percentage in top 25%'); 
text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 11.6], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 11.6])

%end of figure 3

%Figure S4: relating to Figure 3C
load('border conjunctivity shuffled UPDATED hyperfield excluded.mat')
load('border fields in top max firing hyperfield excluded.mat')

figure;
n=3;
m=1;

subplot(n,m,1)
hist(means_20); hold on;
y=0:290;
x= ones(1,291)* mean(border_over_total_20);
plot(x,y,'-','LineWidth', 2, 'Color', [0.85 0.16 0]);
ylabel('shuffles');
xlabel('percentage in top 20%'); 
text(-0.2,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

subplot(n,m,2)
hist(means_30); hold on;
y=0:290;
x= ones(1,291)* mean(border_over_total_30);
plot(x,y,'-','LineWidth', 2, 'Color', [0.85 0.16 0]);
ylabel('shuffles');
xlabel('percentage in top 30%'); 
text(-0.2,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

subplot(n,m,3)
hist(means_40); hold on;
y=0:290;
x= ones(1,291)* mean(border_over_total_40);
plot(x,y,'-','LineWidth', 2, 'Color', [0.85 0.16 0]);
ylabel('shuffles');
xlabel('percentage in top 40%'); 
text(-0.2,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 8.7, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [8.7, 17.4])
%end Figure S4

%Figure 5:Cluster scores
load('cluster score shuffled results.mat')
load('Cluster scores.mat')

figure;
n=1;
m=2;

subplot(n,m,1)
bar(1:length(clusters_bar_graph), clusters_bar_graph); hold on; box off;
ylabel('number of cells');
xlabel('Cluster score'); 
text(-0.2,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

subplot(n,m,2);
hist(mean_cluster_scores_pt_35); hold on; box off;
y= 0:280;
x= ones(1,281)* mean(cluster_scores_35);
plot(x,y,'-','LineWidth', 2);
ylabel('shuffles');
xlabel('mean Cluster score'); 
text(-0.2,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 5.8], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 5.8])

%end Figure 5


% Figure 6

figure;
n=6;
m=5;

cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\stability between fields removed last arena')

load('firing stability rescaled arenas.mat')
rates_corr_rescaled_inds= find(rates_corrcoef >= 0.85);
load('firing stability same arenas.mat')
rates_corr_same_inds= find(rates_corrcoef >= 0.85);
load('firing stability half sessions.mat')
rates_corr_half_inds= find(rates_corrcoef >= 0.85);

% file_names= [{'217_06-07-20.mat_cell5.mat_file3.mat'}; ...
%                 {'214_06-04-21.mat_cell1.mat_file1.mat'};
%                 {'217_06-08-09.mat_cell1.mat_file2.mat'};
%                 {'217_06-08-09.mat_cell1.mat_file1.mat'}];

count=1;



plot_num= [1 2 3 4 7 8 9 10];

letter= [{'A'} {'B'} {'C'} {'D'}];
for h=1:4;

load(file_names{h});
sort_inds=nan(1,length(rates_orig));

    sorted_orig= sort(rates_orig);
    
    for r=1:length(rates_orig)
        sort_inds(r)= find(sorted_orig(r) == rates_orig);
    end
    
    sorted_new=rates_new(sort_inds);
    
    file_name= file_names{h};
    
    subplot(n,m,[plot_num(count) plot_num(count+1)]);
    count=count+2;
    plot(1:length(rates_orig),sorted_orig, 'o-'); hold on;
    plot(1:length(rates_orig),sorted_new, 'ko-'); hold on;
    title(sprintf('%s', file_name(1:end-14)));
   
    text(-0.25,1.05,letter{h},'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
      
end

load('shuffled_remapping_rate_stability 2.mat')
subplot(n,m,[13 14 15]);
hist(rates_corr_mean); hold on;
x= ones(1,301)*0.4323;
y= 0:300;
plot(x,y,'-');
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('shuffled_remapping_rate_stability max field removed2.mat')
subplot(n,m,[16 17 18]);
hist(rates_corr_mean); hold on;
x= ones(1,251)*0.3260;
y= 0:250;
plot(x,y,'-');
text(-0.05,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('fitted arenas stability shuffled.mat')
subplot(n,m,[19 20 21]);
hist(mean_corrcoef); hold on;
x= ones(1,251)*0.361;
y= 0:250;
plot(x,y,'-');
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('fitted arenas stability shuffled wo apex.mat')
subplot(n,m,[22 23 24]);
hist(mean_corrcoef); hold on;
x= ones(1,291)*0.3461;
y= 0:290;
plot(x,y,'-');
text(-0.05,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)


cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS')
load('shuffled_remapping_rate_stability first and last arena2.mat');
subplot(n,m,[25 26 27]);
hist(rates_corr_mean); hold on;
y= 0:270;
x= ones(1,271)*0.499;
plot(x,y,'-');
text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('shuffled_remapping_rate_stability first and last arena wo max2.mat')
subplot(n,m,[28 29 30]);
hist(rates_corr_mean); hold on;
x= ones(1,251)*0.5391;
y= 0:250;
plot(x,y,'-');
text(-0.05,1.05,'F','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

%end Figure 6 

% Figure S3

%end Figure S3

%Figure S5
load('Cluster scores.mat')
load('cluster score shuffled results.mat')

mean_25= nanmean(cluster_scores_25);
mean_45= nanmean(cluster_scores_45); 

figure;
n=1;
m=2;

subplot(n,m,1)
hist(mean_cluster_scores_pt_25); hold on;
y= 0:280;
x= ones(1,281)* mean_25;
plot(x,y,'-','LineWidth', 2, 'Color', [0.85 0.16 0]);
xlabel('mean Cluster score ("near" <= 0.25');
ylabel('shuffles');
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

subplot(n,m,2)
hist(mean_cluster_scores_pt_25); hold on;
y=0:290;
x= ones(1,291)* mean_45;
plot(x,y,'-','LineWidth', 2, 'Color', [0.85 0.16 0]);
xlabel('mean Cluster score ("near" <= 0.45');
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 5.8], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 5.8])
%end Figure S5



