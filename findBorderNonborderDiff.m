function [border_mean_rates, nonborder_means, nonborder_means_wo_max, border_means_wo_max, border_zones] = findBorderNonborderDiff(number_zone_mat, peak_rates_list)


[size_x, size_y] =size(number_zone_mat);

top = unique(number_zone_mat(1,:));
top(find(top==0))= [];

left = unique(number_zone_mat(:,1));
left(find(left==0))= [];

right = unique(number_zone_mat(:,size_y));
right(find(right==0))= [];

bottom = unique(number_zone_mat(size_x,:)) ;
bottom(find(bottom==0))= [];

A = union(top, left);
C = union(bottom, right);

border_zones = union(A, C);
border_zones(find(border_zones==0))= [];

border_mean_rates = peak_rates_list(border_zones);

nonborder_means = setdiff(peak_rates_list, border_mean_rates);

%border_to_nonborder_diff= mean(border_mean_rates)/mean(nonborder_means);

% looks at border to nonborder difference with max peak removed
peak_rate = max(peak_rates_list);

if any(border_mean_rates==peak_rate)
    nonborder_means_wo_max=nonborder_means;
    border_means_wo_max= border_mean_rates;
    border_means_wo_max(border_means_wo_max==peak_rate)=[];
elseif any(nonborder_means ==peak_rate)
    nonborder_means_wo_max=nonborder_means;
    border_means_wo_max= border_mean_rates;
    nonborder_means_wo_max(nonborder_means_wo_max==peak_rate)=[];
else
    disp('ERROR ERROR')
end

%border_to_nonborder_diff_wo_max= mean(border_means_wo_max)/mean(nonborder_means_wo_max);



border_means_wo_max= mean(border_means_wo_max); hold on;
nonborder_means_wo_max= mean(nonborder_means_wo_max); hold on;




%         top_to_central = mean(peak_rates_list(top))/mean(nonborder_means);
%         S.top = top_to_central;
%
%         bottom_to_central = mean(peak_rates_list(bottom))/mean(nonborder_means);
%         S.bottom = bottom_to_central;
%
%         left_to_central= mean(peak_rates_list(left))/mean(nonborder_means);
%         S.left= left_to_central;
%
%         right_to_central= mean(peak_rates_list(right))/mean(nonborder_means);
%         S.right = right_to_central;