function closest_cen_coord= FindClosestCenCoords(max_inds, posi_x, posi_y)

distances= nan(1,length(max_inds));

[max_inds_len,~]=size(max_inds);

for j=1:max_inds_len
    distances(j)= Distance(max_inds(j,1), max_inds(j,2), posi_x, posi_y); 
end

min_ind=find(distances==min(distances));
min_ind=min_ind(1);
closest_cen_coord(1,:)= max_inds(min_ind,:);

disp('')