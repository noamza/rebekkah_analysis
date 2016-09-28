function [median_ordered_distances, peak_score]= SpikeDistanceFromCenter(max_inds, spk_x, spk_y, peak_rates, sorted_means, pos_x, pos_y)

conv_max_inds(:,1)= max_inds(:,1)*5; %- max(pos_x);
conv_max_inds(:,2)= max_inds(:,2)*5 ;%- max(pos_y);
    
distances_from_peaks=[];
    
distances_spks_from_peak=nan(length(max_inds), length(spk_x)); 


    for h=1:length(spk_x)
        for j= 1:length(max_inds)
            distances_from_peaks(j)= Distance(spk_x(h), spk_y(h), conv_max_inds(j,2), conv_max_inds(j,1)); %find distance of spk from zones
        end
        min_dist_zone_num= find(distances_from_peaks==min(distances_from_peaks)); %finds which zone spk belongs to
        distances_spks_from_peak(min_dist_zone_num, h)= min(distances_from_peaks);    
end
    
    for h=1:length(peak_rates)
        max_zone_num(h)= find((peak_rates==sorted_means(h))); %order of number zones from lowest to highest firing rate
    end
    
    mean_distances_spks_from_peaks= nanmean(distances_spks_from_peak.');
    median_distances_spks_from_peaks= nanmedian(distances_spks_from_peak.');
    std_distances_spks_from_peaks= nanstd(distances_spks_from_peak.');
    
    %order mean distances by lowest t strongest firing rate fields
    
    for h= 1:length(max_zone_num)
        ordered_distances(h)= mean_distances_spks_from_peaks(max_zone_num(h));
        median_ordered_distances(h)= median_distances_spks_from_peaks(max_zone_num(h));
        std_ordered_distances(h)= std_distances_spks_from_peaks(max_zone_num(h));
    end
    
    
    %peak score should all be above 1 
    peak_score= sorted_means(end)/mean(sorted_means(1:end-1)); % how much higher the peak is from the mean of all the fields;
    
    





