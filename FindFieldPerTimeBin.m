function field_t = FindFieldPerTimeBin(pos_t, pos_x_inds, pos_y_inds, number_zone_mat) 


for time=1:length(pos_t)
    if isnan(pos_x_inds(time))
        field_t(time)=nan;
    else
        field_t(time)= number_zone_mat(pos_x_inds(time), pos_y_inds(time)); %number of zonemat at each timestamp
    end
end