cd('N:\users\rebekkah\results and info of analysis');

figure;
n=5;
m=6;

load('firing stability half sessions.mat')
rates_corr_half= rates_corrcoef;
rates_corr_half_inds= find(rates_corrcoef >= 0.9);
rates_corr_half_inds(2)= [];
load('firing stability rescaled arenas.mat')
rates_corr_rescaled= rates_corrcoef;
rates_corr_rescaled_inds= find(rates_corrcoef >= 0.95);
rates_corr_rescaled_inds(1:5)= [];
rates_corr_rescaled_inds(2:3)= [];
load('firing stability same arenas.mat')
rates_corr_same= rates_corrcoef;
rates_corr_same_inds= find(rates_corrcoef >= 0.85);
rates_corr_same_inds(1:7)= [];


count=1;

plot_num= [1 2 3 4 7 8 9 10];

letter= [{'A'} {'B'} {'C'} {'D'}];

inds= rates_corr_half_inds;

a=1;
b=1;
for h=1:4;
        
    if h==1 || h==2 
        cell_num= rates_corr_half_inds(h);
        
        rates_corr_half(cell_num);
        
        load('peak rates of half sessions.mat')
        load('file names half sessions.mat')
        
        peaks_one= peak_rates_b{cell_num};
        peaks_two= peak_rates_e{cell_num};
        
          file_name= file_names{cell_num};
          file_name=file_name(9:end-4);
        
    elseif h==3 || h==4 
        cell_num= rates_corr_rescaled_inds(b);
        
        
        rates_corr_rescaled(cell_num);
        
        load('firing stability rescaled arenas.mat')
        
        
        peaks_one= rates_orig{cell_num};
        peaks_two= rates_next{cell_num};
        
        load('filenames of rescaling firing comparison.mat')
        file_name= filename_rescaling{cell_num};
        
        b=b+1;
    end
    
    
    %remove all zero and nans from both rate vectors
    zero_inds_1= find(peaks_one==0);
    zero_inds_2= find(peaks_two==0);
    nan_inds_1= find(isnan(peaks_one));
    nan_inds_2= find(isnan(peaks_two));
    zero_inds= union(zero_inds_1, zero_inds_2);
    nan_inds= union(nan_inds_1, nan_inds_2);
    remove_inds= union(zero_inds,nan_inds);
    peaks_one(remove_inds)=[];
    peaks_two(remove_inds)=[];
    
    sort_inds=nan(1,length(peaks_one));
    sorted_orig= sort(peaks_one);
    
    for r=1:length(peaks_one)
        sort_inds(r)= find(sorted_orig(r) == peaks_one);
    end
    
    sorted_new=peaks_two(sort_inds);
    
    
    subplot(n,m,[plot_num(count) plot_num(count+1)]);
    count=count+2;
    plot(1:length(peaks_one),sorted_orig, 'o-'); hold on;
    plot(1:length(peaks_one),sorted_new, 'ko-'); hold on;
    title(sprintf('%s', file_name), 'Interpreter', 'none');
    box off;
    
    text(-0.25,1.05,letter{h},'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    
end

load('firing stability half sessions shuffled.mat')
load('firing stability half sessions.mat')
subplot(n,m,[13 14 15]);
hist(mean_corrcoef_half_sessions); hold on;
x= ones(1,301)* mean(rates_corrcoef);
y= 0:300;
plot(x,y,'-');
text(-0.15,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

subplot(n,m,[16 17 18]);
hist(mean_corrcoef_half_sessions_wo_max); hold on;
x= ones(1,301)*mean(rates_corrcoef_wo_max);
y= 0:300;
plot(x,y,'-');
text(-0.15,1.05,'F','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

load('firing stability rescaling data shuffled.mat')
load('firing stability rescaled arenas.mat')
subplot(n,m,[19 20 21]);
hist(mean_corrcoef_rescaled_arenas); hold on;
x= ones(1,301)* mean(rates_corrcoef);
y= 0:300;
plot(x,y,'-');
text(-0.15,1.05,'G','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

subplot(n,m,[22 23 24]);
hist(mean_corrcoef_rescaled_arenas_wo_max); hold on;
x= ones(1,301)*nanmean(rates_corrcoef_wo_max);
y= 0:300;
plot(x,y,'-');
text(-0.15,1.05,'H','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

load('firing stability same arenas.mat')
subplot(n,m,[25 26 27]);
hist(mean_corrcoef_same_arenas); hold on;
y= 0:300;
x= ones(1,301)*mean(rates_corrcoef);
plot(x,y,'-');
text(-0.15,1.05,'I','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

subplot(n,m,[28 29 30]);
hist(mean_corrcoef_same_arenas_wo_max); hold on;
x= ones(1,301)*mean(rates_corrcoef_wo_max);
y= 0:300;
plot(x,y,'-');
text(-0.15,1.05,'G','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
box off;

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 29], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 29])
