
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat',...
    'number_zone_mat', 'max_indices_ns', 'peak_rates_ns_all')

border_mean= nan(1,length(peak_rates_ns_all));
central_mean= nan(1,length(peak_rates_ns_all));

for i=1:length(peak_rates_ns_all)
 
    peak_rates=peak_rates_ns_all{i};
    max_inds=max_indices_ns{i};
    
    hyper_ind=find(peak_rates==max(peak_rates));
    peak_rates(hyper_ind)=[];
    max_inds(hyper_ind,:)=[];
    
[border_mean_rates, nonborder_means] = BorderNonborderDiff(number_zone_mat{i}, max_inds, peak_rates);
    
border_mean(i)=  mean(border_mean_rates);
central_mean(i)= mean(nonborder_means);

end

wilcoxon_p_value= signrank(border_mean,central_mean);

save('border vs center mean rates ns wo hyper', 'border_mean', 'central_mean', 'wilcoxon_p_value')

%plots results
figure; 
x_lim=max([max(border_mean) max(central_mean)]);
plot(border_mean, central_mean, 'ko'); hold on;
axis equal;
x=0:x_lim+1;
plot(x,x,'k-');
xlim([0 x_lim])

figure;
x=border_mean;
stderr_b=nanstd(x)/sqrt(length(x));
x=central_mean;
stderr_c=nanstd(x)/sqrt(length(x));
errorbar([nanmean(border_mean) nanmean(central_mean)],[stderr_b stderr_c])