function [dist] = FindDistTwoHypers(hyper_ind_1, hyper_ind_2,...
    max_inds_2)

% find two hyperfield dist diff
dist= Distance(hyper_ind_1(1),hyper_ind_1(2),hyper_ind_2(1),hyper_ind_2(2));

[lenn,~]= size(max_inds_2);

% find min dist
field_dists=nan(1,lenn);
for k=1:lenn
    field_dists(k)= Distance(hyper_ind_1(1),hyper_ind_1(2),...
        max_inds_2(k,1), max_inds_2(k,2));
end

min_dist=min(field_dists);

dist= dist- min_dist; 
