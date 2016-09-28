function [border_mean_rates, nonborder_means, border_zones] = BorderNonborderDiff(number_zone_mat, max_inds, peak_rates_list)

[size_x, size_y]= size(number_zone_mat);

for h= 1:length(max_inds)
    [~, norm_dist]= findDistPtToBorder([1 size_x], [1 size_y], max_inds(h,:));
    
    if norm_dist<=0.1
        location(h)= 'b';
    elseif norm_dist >0.1
        location(h)='c';
    end
end

border_zones= find(location=='b');

%to remove max field: add
%max_inds=find(max(peak_rates));


border_mean_rates = peak_rates_list(border_zones);

nonborder_means = setdiff(peak_rates_list, border_mean_rates);

