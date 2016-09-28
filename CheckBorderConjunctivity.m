function [border_over_total]= CheckBorderConjunctivity(number_of_PF, sorted_means, number_zone_mat, peak_rates_list, max_inds,percentage) 

check_top= round(number_of_PF* percentage); % how much is top% of field number (rounded)

[~,~,border_zones] = BorderNonborderDiff(number_zone_mat, max_inds, peak_rates_list); %get border inds with peak_rates (un-sorted)

sorted_ind= nan(1,length(peak_rates_list));
for h=1:length(peak_rates_list)
    ind=find(sorted_means(h)==peak_rates_list);
    if length(ind)>1            %if there are two of the same values
        ind=ind(1);             % choose first
        peak_rates_list(ind)=peak_rates_list(ind)+0.01; %change it so next run it will catch second same value 
    end
    sorted_ind(h)= ind;
end

sorted_border(1:length(peak_rates_list))= 'c';
sorted_border(border_zones)= 'b';
sorted_border= sorted_border(sorted_ind); %order 'b' or 'c' in increasing order by firing rate

sum_border= sum(sorted_border(end-check_top+1:end)=='b');
sum_total= length(sorted_border(end-check_top+1:end));

border_over_total= sum_border/ sum_total;

disp('');