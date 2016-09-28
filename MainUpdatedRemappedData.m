function [fano_factor, norm_max_index, max_indices, ...
    max_indices_sm, norm_dist, rm_size, peak_rates_all, ...
    peak_rates_all_sm, norm_second_dist, number_zone_mat_sm,...
    number_zone_mat]= findAllInfo(pos_mean_x,pos_mean_y,pos_t,spk_x,spk_y,spk_t)

parms.bin_size=6;
rate_mat_ns= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,pos_t,spk_x,spk_y,spk_t,parms);

% Create AutoCorr
autocorr=Cross_Correlation(rate_mat, rate_mat);
auto_max_inds= FindAutoMaxInds(autocorr);
PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);

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
norm_max_index= max_index ./ size_rate_mat;

max_indices= max_inds_ns; %save for future max peak shuffling
max_indices_sm= max_inds; %save for future max peak shuffling

%find distance of max field location to border
[~,norm_dist]= findDistPtToBorder([1 size_rate_mat(1)], [1 size_rate_mat(2)], max_index);

%save rate mat sizes and nonsmooth peak rates for future max peak shuffling
rm_size= size_rate_mat;
peak_rates_all=peak_rates_ns;

peak_rates_all_sm=peak_rates; %save for border vs nonborder mean rates analysis

size_rate_mat= size(rate_mat);
%     %find distance of second max field location to border
peak_rates_wo_max= peak_rates;
peak_rates_wo_max(peak_rates_wo_max==max(peak_rates_wo_max))= nan;
second_ind= find(peak_rates_wo_max==max(peak_rates_wo_max));
second_ind=second_ind(1); %in case same rate take first
second_max_index= max_inds(second_ind,:);
[~,norm_second_dist]= findDistPtToBorder([ 1 size_rate_mat(1)], [1 size_rate_mat(2)], second_max_index);

% create zone mats for images and number_zone_mat for future shuffling
[~, number_zone_mat_sm]= CreateZoneMat(rate_mat, PF_radius, max_inds, peak_rates);
[~, number_zone_mat]= CreateZoneMat(rate_mat_ns, PF_radius_ns, max_inds_ns, peak_rates_ns);

%image
%     fig=figure;
%     n=2;m=4;
%     subplot(n,m,1)
%     plot(pos_mean_x,pos_mean_y,'k');hold on;
%     plot(spk_x,spk_y,'.r');
%     axis equal;axis off;
%     axis ij;
%     title_name= file_name(14:end-4);
%     title(sprintf('%s', title_name), 'Interpreter', 'none');
%
%     subplot(n,m,2)
%     imagesc(rate_mat);
%     axis equal; axis off;
%     title(sprintf('%0.1f Hz', max(peak_rates)), 'HorizontalAlignment', 'left');
%
%     subplot(n,m,3)
%     imagesc(zone_mat);
%     axis equal; axis off;
%
%     subplot(n,m,4)
%     plot(1:length(peak_rates),sort(peak_rates), 'ko-');
%     title(sprintf('Fano factor=%0.1f', fano_factor(i)));
%
%     subplot(n,m,5)
%     imagesc(autocorr);
%     axis equal; axis off;
%
%     subplot(n,m,6)
%     imagesc(rate_mat_ns);
%     axis equal; axis off;
%     title(sprintf('%0.1f Hz', max(peak_rates_ns)),'HorizontalAlignment', 'left');
%
%     subplot(n,m,7)
%     imagesc(zone_mat_ns);
%     axis equal; axis off;
%     title(sprintf('dist=%0.2f', norm_dist(i)));
%
%     subplot(n,m,8)
%     plot(1:length(peak_rates_ns),sort(peak_rates_ns), 'ko-');
%     title(sprintf('Fano factor=%0.1f', fano_factor_ns));
%
%     %save images
%     cd(parms.dir_save_pictures);
%     saveas(fig,sprintf('%s.jpg',file_name)); %
%     cd(parms.dir_load_data);

close all

disp('');
