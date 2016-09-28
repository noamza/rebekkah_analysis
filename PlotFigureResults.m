
fig=figure;

n=5;
m=6;

%3A: 2nd max peak
%subplot(n,m,1)
%cd('N:\users\rebekkah\bin size 6 nonsmooth');
%cd('N:\users\rebekkah\final data smoothed')
%cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\stability between fields removed last arena')
% 
% 
% 
% file_names= [{'217_06-07-20.mat_cell5.mat_file3.mat'}; ...
%                 {'214_06-04-21.mat_cell1.mat_file1.mat'};
%                 {'217_06-08-09.mat_cell1.mat_file2.mat'};
%                 {'217_06-08-09.mat_cell1.mat_file1.mat'}];
%count=1;

%plot_num= [1 2 3 4 7 8 9 10];


for h=1:4;

load(file_names{h});

sorted_orig=[];
sort_inds=[];
sorted_new=[];

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
%   
    if h==1
        text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    end
    
   
end

% figure;
n=2;
m=1;
cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS')


%%% figure 2A-E
% figure;
% n=3;
% m=2;
% 
% cd('N:\users\rebekkah\final data smoothed');
% load('fano factors.mat');
% subplot(n,m,1);
% hist(fano_factors); hold on;
% text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% cd('N:\users\rebekkah');
% load('simulated results.mat');
% subplot(n,m,3);
% % fix bug later
% fanos_med(fanos_med>1)= [];
% hist(fanos_med); hold on;
% y= 0:200;
% x= ones(1,201)* 2.186; 
% plot(x,y,'-');
% text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% cd('N:\users\rebekkah\bin size 6 nonsmooth');
% load('distances PF 8 bin6 nonsmooth.mat')
% subplot(n,m,2)
% hist(distances); hold on;
% text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% load('shuffled percentages.mat')
% subplot(n,m,4)
% hist(max_sh); hold on;
% y= 0:300;
% x= ones(1,301)* (56/86); 
% plot(x,y,'-');
% text(-0.25,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% cd('N:\users\rebekkah');
% load('simulated results.mat');
% subplot(n,m,6);
% hist(sums_ns); hold on;
% y= 0:200;
% x= ones(1,201)* 0.6512; 
% plot(x,y,'-');
% text(-0.25,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)


% load('shuffled_remapping_rate_stability 2.mat')
% subplot(n,m,[13 14 15]);
% hist(rates_corr_mean); hold on;
% x= ones(1,301)*0.4323;
% y= 0:300;
% plot(x,y,'-');
% text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% load('shuffled_remapping_rate_stability max field removed2.mat')
% subplot(n,m,[16 17 18]);
% hist(rates_corr_mean); hold on;
% x= ones(1,251)*0.3260;
% y= 0:250;
% plot(x,y,'-');
% %text(-0.05,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% load('fitted arenas stability shuffled.mat')
% subplot(n,m,[19 20 21]);
% hist(mean_corrcoef); hold on;
% x= ones(1,251)*0.361;
% y= 0:250;
% plot(x,y,'-');
% text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% load('fitted arenas stability shuffled wo apex.mat')
% subplot(n,m,[22 23 24]);
% hist(mean_corrcoef); hold on;
% x= ones(1,291)*0.3461;
% y= 0:290;
% plot(x,y,'-');
% %text(-0.05,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% 
% cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS')
% load('shuffled_remapping_rate_stability first and last arena2.mat');
% subplot(n,m,[25 26 27]);
% hist(rates_corr_mean); hold on;
% y= 0:270;
% x= ones(1,271)*0.499;
% plot(x,y,'-');
% text(-0.25,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% load('shuffled_remapping_rate_stability first and last arena wo max2.mat')
% subplot(n,m,[28 29 30]);
% hist(rates_corr_mean); hold on;
% x= ones(1,251)*0.5391;
% y= 0:250;
% plot(x,y,'-');
% %text(-0.05,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% set(fig, 'units', 'centimeters', 'Position', [0 0 17.4 22.8])
% 
% % subplot(n,m,[14 15 19 20]);
% % hist(rates_corr_sum2); hold on;
% % x= ones(1,271)*27;
% % y= 0:270;
% % plot(x,y,'-');
% % text(-0.05,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% %  
% %     load('Results_Cell_r11138_d200405_s02_t6_c1.mat')
% %     
% %     pos_mean_x=(S.pos.x + S.pos.x2)/2;
% %     pos_mean_y=(S.pos.y + S.pos.y2)/2;
% %     
% %     % build the axis of location when spikes are made
% %     spk_x=interp1(S.pos.t,pos_mean_x,S.spk.t);
% %     spk_y=interp1(S.pos.t,pos_mean_y,S.spk.t);
% %     
% %     plot(pos_mean_x,pos_mean_y,'k');hold on;
% %     plot(spk_x,spk_y,'.r');
% %     axis square; axis off;
% %     axis ij
% %     title('Cell_r11138_d200405_s02_t6_c1');
% %     
% %     hold on;
%     
%     
%     
%     %hist(cluster_scores); hold on;
%     
%     %load('second max peak distances from border stats.mat');
%     %load('distances PF 8 bin6 nonsmooth.mat')
%     
%     % cut_off= 0.05;
%     % y_lim_mark= 270;
%     % real_value= sum(distances<= cut_off);
%     % x=ones(1,y_lim_mark+1) * real_value;
%     % y=0:y_lim_mark;
%     % plot(x,y,'--', 'Linewidth', 2.5, 'color', 'r');
%     %
%     text(-0.2,1.05,'A)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
%     
%     %................
%     
%     subplot(n,m,2)
%     %load('mean cluster score 1000x.mat')
%     
%     %hist(total_mean_cluster); hold on;
%     
%     imagesc(S.rate_mat); hold on;
%     axis off; axis square;
%     title(sprintf('%0.2f', max(max(S.rate_mat))));
%     
%     %cut_off= 0.15;
%     % y_lim_mark= 270;
%     % real_value= mean(cluster_scores);
%     % x=ones(1,y_lim_mark+1) * real_value;
%     % y=0:y_lim_mark;
%     % plot(x,y,'--', 'Linewidth', 2.5, 'color', 'r');
%     
%     text(-0.2,1.05,'B)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
%     
%     %3C: top quarter border amount
%     % cd('N:\users\rebekkah\final data smoothed');
%     % subplot(n,m,3)
%     % load('border top 25 percent 1000x.mat');
%     % hist(means); hold on;
%     %
%     % load('border over total 25 percent.mat');
%     % cut_off= mean(border_over_total);
%     % y_lim_mark= 270;
%     % real_value= cut_off;
%     % x=ones(1,y_lim_mark+1) * real_value;
%     % y=0:y_lim_mark;
%     % plot(x,y,'--', 'Linewidth', 3, 'color', 'r');
%     %
%     % text(-0.3,0.98,'C)','Units', 'Normalized', 'VerticalAlignment', 'Top')
%     %
%     % %3B: border vs central scatterplot
%     % subplot(n,m,2)
%     % load('border vs center mean rates.mat');
%     % scatter(border_mean, central_mean); hold on;
%     % axis equal;
%     %
%     % x= 0:40;
%     % y=0:40;
%     % plot(x,y, '--', 'color', 'r');
%     %
%     % text(-0.3,0.98,'B)','Units', 'Normalized', 'VerticalAlignment', 'Top')
%     
%     subplot(n,m,3);
%     imagesc(S.zone_mat);
%     axis off; axis square;
%     text(-0.2,1.05,'C)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw', 'b', 'fontsize', 14)
%     
%     subplot(n,m,4);
%     imagesc(S.autocorr);
%     axis off; axis square;
%     text(-0.2,1.05,'D)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw', 'b', 'fontsize', 14)
%     
%     subplot(n,m,5);
%     h= 1:length(S.sorted_means);
%     plot(h, S.sorted_means);
%     axis square;
%     text(-0.2,1.05,'E)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw', 'b', 'fontsize', 14)


%% FIGURE 6

cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\stability between fields removed last arena') 

file_names= [{'217_06-07-20.mat_cell5.mat_file3.mat'}; ...
                 {'214_06-04-21.mat_cell1.mat_file1.mat'};
                 {'217_06-08-09.mat_cell1.mat_file2.mat'};
                 {'217_06-08-09.mat_cell1.mat_file1.mat'}];
count=1;

plot_num= [1 2 3 4 7 8 9 10];

for h=1:4;

load(file_names{h});

sorted_orig=[];
sort_inds=[];
sorted_new=[];

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
%   
    %if h==1
        text(-0.05,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    %end
    
   
end

cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS')

load('shuffled_remapping_rate_stability 2.mat')
subplot(n,m,[13 14 15]);
hist(rates_corr_mean); hold on;
x= ones(1,301)*0.4323;
y= 0:300;
plot(x,y,'-');
text(-0.05,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('shuffled_remapping_rate_stability max field removed2.mat')
subplot(n,m,[16 17 18]);
hist(rates_corr_mean); hold on;
x= ones(1,251)*0.3260;
y= 0:250;
plot(x,y,'-');
text(-0.05,1.05,'H','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('fitted arenas stability shuffled.mat')
subplot(n,m,[19 20 21]);
hist(mean_corrcoef); hold on;
x= ones(1,251)*0.361;
y= 0:250;
plot(x,y,'-');
text(-0.05,1.05,'F','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('fitted arenas stability shuffled wo apex.mat')
subplot(n,m,[22 23 24]);
hist(mean_corrcoef); hold on;
x= ones(1,291)*0.3461;
y= 0:290;
plot(x,y,'-');
text(-0.05,1.05,'I','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

cd('N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS')
load('shuffled_remapping_rate_stability first and last arena2.mat');
subplot(n,m,[25 26 27]);
hist(rates_corr_mean); hold on;
y= 0:270;
x= ones(1,271)*0.499;
plot(x,y,'-');
text(-0.05,1.05,'G','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

load('shuffled_remapping_rate_stability first and last arena wo max2.mat')
subplot(n,m,[28 29 30]);
hist(rates_corr_mean); hold on;
x= ones(1,251)*0.5391;
y= 0:250;
plot(x,y,'-');
text(-0.05,1.05,'J','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

set(fig, 'units', 'centimeters', 'Position', [0 0 17.4 22.8])

% subplot(n,m,[14 15 19 20]);
% hist(rates_corr_sum2); hold on;
% x= ones(1,271)*27;
% y= 0:270;
% plot(x,y,'-');
% text(-0.05,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)



set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 14], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 14])