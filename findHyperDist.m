function [norm_dist, peak_rates_all]= findHyperDist(rate_mat,PF_radius)

PF_radius_ns= PF_radius/2; %half since bin_size is 6 instead of 3
max_inds_ns= FindMaxIndsRateMap(rate_mat);
max_inds_ns= RemoveTooCloseMaxInds(max_inds_ns, PF_radius_ns, rate_mat, 1.7);

[len,~]=size(max_inds_ns);

% Find peak firing rate at max_inds
peak_rates_ns=findPeakRates(max_inds_ns, rate_mat); 

if len>=3
% Find location of max-firing field (hyperfield)- nonsmoothed data
ind= find(peak_rates_ns==max(peak_rates_ns));
ind=ind(1); %in cases where same rate, take first
max_index= max_inds_ns(ind,:);
size_rate_mat= size(rate_mat);


%find distance of max field location to border
[~,norm_dist]= findDistPtToBorder([1 size_rate_mat(1)], [1 size_rate_mat(2)], max_index);

else 
    norm_dist=nan;
end

%save rate mat sizes and nonsmooth peak rates for future max peak shuffling
peak_rates_all=peak_rates_ns;

disp('');
