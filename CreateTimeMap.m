function [time_mat]=CreateTimeMap(posx,posy,post,parms)
%function [rate_mat, spike_mat_smooth]=CreateRateMap(posx,posy,post,spkx,spky,spkt,parms)
max_x = (max(posx)); 
max_y = (max(posy));
min_x = min((posx));
min_y = min((posy));

% divid the environment into spatial bins 
axis_x = min_x:parms.bin_size:max_x;
axis_y = min_y:parms.bin_size:max_y;
dt=post(2)-post(1);

time_mat = zeros(length(axis_y),length(axis_x));
spike_mat = zeros(length(axis_y),length(axis_x));
rate_mat = zeros(length(axis_y),length(axis_x));

%create time mat (compute how much time the rat spends in each bin)
% find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to
% 
for i = 1:length(post)
    if ~isnan(posx(i)) && ~isnan(posy(i))
        [~,x_ind] =  min(abs(posx(i)-axis_x));
        [~,y_ind] =  min(abs(posy(i)-axis_y));
        time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind)+dt;
        
    end
end


time_mat=SmoothRateMat(time_mat,parms);
 
% spike_mat_smooth =SmoothRateMat(spike_mat, parms);  % want to see spike_mat irrelevant to time 
 
 disp('');
 