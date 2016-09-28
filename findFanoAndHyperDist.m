function [fano_factor, norm_dist, peak_rates_all_sm]= findFanoAndHyperDist(pos_mean_x,pos_mean_y,pos_t,spk_x,spk_y,spk_t,rate_mat,PF_radius)

parms.bin_size=6;
rate_mat_ns= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,pos_t,spk_x,spk_y,spk_t,parms);

% Find Max_Inds of smoothed & nonsmoothed rate mat
max_inds= FindMaxIndsRateMap(rate_mat);
strength= 1.9;
max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength);

PF_radius_ns= PF_radius/2; %half since bin_size is 6 instead of 3
max_inds_ns= FindMaxIndsRateMap(rate_mat_ns);
max_inds_ns= RemoveTooCloseMaxInds(max_inds_ns, PF_radius_ns, rate_mat_ns, strength);

% Find peak firing rate at max_inds
peak_rates= nan(1,length(max_inds));
for cen= 1:length(max_inds);
    peak_rates(cen)= rate_mat(max_inds(cen,1), max_inds(cen,2));
end

peak_rates_ns= nan(1,length(max_inds_ns));
for cen= 1:length(max_inds_ns);
    peak_rates_ns(cen)= rate_mat_ns(max_inds_ns(cen,1), max_inds_ns(cen,2));
end

% Find Fano factor
fano_factor= var(peak_rates)/mean(peak_rates);
%     fano_factor_ns= var(peak_rates_ns)/mean(peak_rates);
%
% Find location of max-firing field (hyperfield)- nonsmoothed data
ind= find(peak_rates_ns==max(peak_rates_ns));
ind=ind(1); %in cases where same rate, take first
max_index= max_inds_ns(ind,:);
size_rate_mat= size(rate_mat_ns);

%find distance of max field location to border
[~,norm_dist]= findDistPtToBorder([1 size_rate_mat(1)], [1 size_rate_mat(2)], max_index);

%save rate mat sizes and nonsmooth peak rates for future max peak shuffling
peak_rates_all=peak_rates_ns;

peak_rates_all_sm=peak_rates; %save for border vs nonborder mean rates analysis

close all

disp('');
