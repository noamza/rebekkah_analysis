

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results sorted PF 93 cells';
%parms.dir_load_data2 = 'N:\users\rebekkah\final data smoothed\data sets\results spikemat 202 cells';
%parms.dir_save_pictures = 'N:\users\rebekkah\final data smoothed\data sets\3 spike and rate mats images 93 cells';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for shuffle_x_times=1:500

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};
    dat = load(file_name);
    RateResults= dat.S;
    
%     cd(parms.dir_load_data2);
%     dat2 = load(file_name);
%     SmoothSpike= dat2.S;
    
    %finds the percentage of PF size of maximal firing field
    max_field_rate= RateResults.sorted_means(end);
    PF_areas= RateResults.PF_areas(randperm(length(RateResults.PF_areas)));
    max_rate_field_size=PF_areas(end);
    percentage_field(k) = max_rate_field_size/max(RateResults.PF_areas);
    
    %finds the percentage of PF size of maximal firing field once
    %multiplied by the percentage of its size
    
%     [size_accounted_max_field_rate, index] = max(RateResults.rate_over_area);
%     size_accounted_max_rate_field_size = PF_areas(index);
%     size_accounted_percentage_field(k)=size_accounted_max_rate_field_size/max(RateResults.PF_areas);
%     
%     %finds the number of spikes per field (instead of the firing rate)
%     spike_rates=unique(RateResults.zone_spike_mat);
%     spike_rates(spike_rates==0)= [];
%     
%     %sorts the spikes
%     sorted_spike_rates= sort(spike_rates)';
%     
%     max_spike_rate= sorted_spike_rates(end);
%     
% %     spike_PF_areas=[];
%     
%     for h=1:length(sorted_spike_rates)
%         spike_PF_area= sum(RateResults.zone_spike_mat==sorted_spike_rates(h));
%         spike_PF_area=sum(spike_PF_area);
%         spike_PF_areas(h)= spike_PF_area; % areas of place fields in same order as sorted means
%     end
%     
%     spike_PF_areas= spike_PF_areas(randperm(length(spike_PF_areas)));
%     
%     max_spike_field_size=spike_PF_areas(end);
%     percentage_field_spikes(k) = max_spike_field_size/max(spike_PF_areas);
%     
%     max_PF_spike_area= max(spike_PF_areas); %maximum PF size to be used at standard PF size
%     rate_over_area_spike= sorted_spike_rates./(spike_PF_areas/max_PF_spike_area); %divide number by percentage of PF size
%     
%     [size_accounted_max_field_spike, Index] = max(rate_over_area_spike);
%     size_accounted_max_rate_field_spike = spike_PF_areas(Index);
%     size_accounted_percentage_field_spikes(k)=size_accounted_max_rate_field_spike/max(spike_PF_areas);
%     
%     
%     %finds the number of spikes per field using pure spike mat (not
%     %derived from rate mat)
%     
%     pure_spike_rates=unique(SmoothSpike.zone_mat);
%     pure_spike_rates(pure_spike_rates==0)= [];
%     
%     %sorts the spikes
%     pure_sorted_spike_rates= sort(pure_spike_rates)';
%     
%     pure_max_spike_rate= pure_sorted_spike_rates(end);
%     
%     pure_spike_PF_areas=[];
%     
%     for h=1:length(pure_sorted_spike_rates)
%         pure_spike_PF_area= sum(SmoothSpike.zone_mat==pure_sorted_spike_rates(h));
%         pure_spike_PF_area=sum(pure_spike_PF_area);
%         pure_spike_PF_areas(h)= pure_spike_PF_area; % areas of place fields in same order as sorted means
%     end
%     
%     pure_spike_PF_areas= pure_spike_PF_areas(randperm(length(pure_spike_PF_areas)));
%     
%     pure_max_spike_field_size=pure_spike_PF_areas(end);
%     pure_percentage_field(k) = pure_max_spike_field_size/max(pure_spike_PF_areas);
%     
%     pure_max_PF_spike_area= max(pure_spike_PF_areas); %maximum PF size to be used at standard PF size
%     pure_rate_over_area= pure_sorted_spike_rates./(pure_spike_PF_areas/pure_max_PF_spike_area); %divide number by percentage of PF size
%     
%     [pure_size_accounted_max_field, Index] = max(pure_rate_over_area);
%     pure_size_accounted_max_rate_field = pure_spike_PF_areas(Index);
%     pure_size_accounted_percentage_field(k)=pure_size_accounted_max_rate_field/max(pure_spike_PF_areas);
    
%     n=1; m=3;
%     fig=figure; 
%     subplot(n,m,1)
%     imagesc(RateResults.zone_mat)
% 
%      subplot(n,m,2)
%     imagesc(RateResults.zone_spike_mat)
%      subplot(n,m,3)
%     imagesc(SmoothSpike.zone_mat)
  
   % cd(parms.dir_save_pictures);
    %  saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
   % saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',RateResults.rat,RateResults.date,RateResults.session,RateResults.tetrode,RateResults.cell)); %
    
end

percentage_field_sums(shuffle_x_times)=sum(percentage_field>0.5 & percentage_field< 0.6); 

percentage_field_below_point_6(shuffle_x_times)=sum(percentage_field< 0.6); 

end

save('percentage_field', 'percentage_field_sums');
save('percentage_field_below_point_6', 'percentage_field_below_point_6');

figure;hist(percentage_field_sums);
figure;hist(percentage_field_below_point_6);