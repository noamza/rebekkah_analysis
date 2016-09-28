%pos_time is timestamp at each timebin in vector (Cell.pos.t)
%pos_x_inds and pos_y_inds are coordinates of trajectory converted to rate
%map coordinates

function firing_rate_t = FindFiringRatePerTimeBin(pos_time, pos_x_inds, pos_y_inds, rate_mat) 

%creates vector of firing rate at each time bin 
firing_rate_t=nan(1,length(pos_time)); 
for time=1:length(pos_time)
            if isnan(pos_x_inds(time))
                firing_rate_t(time)=nan;
            else
                firing_rate_t(time)= rate_mat(pos_x_inds(time), pos_y_inds(time)); %expected firing rate at each timestamp derived using ratemat
            end
end

%firing_rate_t is in units of spikes/sec (same at rate_mat)