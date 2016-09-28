function [six_orientation_pts] = find_six_points(autocorr)

%find autocorrelation mat maximum points
[auto_max_inds] = FindMaxIndsRateMap(autocorr);

%finds distances from center of auto_corr_map to all max peaks
cen = 1:length(auto_max_inds);
auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);

%find point closest to center
min_distance_index = find(auto_distances==min(auto_distances));
middle_pt(1) = auto_max_inds(min_distance_index,1);
middle_pt(2) = auto_max_inds(min_distance_index,2);

%repeats find distances with more accurate center point
cen = 1:length(auto_max_inds);
auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),middle_pt(1),middle_pt(2));

[new_distances, six_inds] = sort(auto_distances (:));
%% finds the points for the module ellipse
% finds the 6 closest pts to the center
% auto_dist_inds1 = find(auto_distances == new_distances(2));
% auto_dist_inds2 = find(auto_distances == new_distances(3));
% auto_dist_inds3 = find(auto_distances == new_distances(4));
% auto_dist_inds4 = find(auto_distances == new_distances(5));
% auto_dist_inds5 = find(auto_distances == new_distances(6));
% auto_dist_inds6 = find(auto_distances == new_distances(7));
% 
% union1 = union(auto_dist_inds1, auto_dist_inds2);
% union2 = union(auto_dist_inds3, auto_dist_inds4);
% union3 = union(auto_dist_inds5, auto_dist_inds6);
% union1 = union(union1, union2);
% 
% auto_dist_inds = union(union1, union3);

six_orientation_pts=nan(length(auto_dist_inds),2);
for k= 1:length(auto_dist_inds);
    six_orientation_pts (k,1) = auto_max_inds(auto_dist_inds(k),1);
    six_orientation_pts (k,2) = auto_max_inds(auto_dist_inds(k),2);
end

%adds center point

six_orientation_pts (7,1) = middle_pt(1);
six_orientation_pts (7,2) = middle_pt(2);

% to check for accuracy:
% figure; imagesc(autocorr); hold on;
% plot(six_orientation_pts(:,2), six_orientation_pts(:,1), 'x')

disp('');
