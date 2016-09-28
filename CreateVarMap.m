function [var_map]=CreateVarMap(posx,posy,post, parms,firing_var_t, rate_mat)
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
var_sum_mat = zeros(length(axis_y),length(axis_x));

%create time mat (compute how much time the rat spends in each bin)
% find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to

for i = 1:length(post)
    if ~isnan(posx(i)) && ~isnan(posy(i))
        [~,x_ind] =  min(abs(posx(i)-axis_x));
        [~,y_ind] =  min(abs(posy(i)-axis_y));
        time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind)+dt;
    end
end

%create variability map
for i = 1:length(post)
    if ~isnan(posx(i)) && ~isnan(posy(i))
        [~,x_ind] =  min(abs(posx(i)-axis_x));
        [~,y_ind] =  min(abs(posy(i)-axis_y));
        var_sum_mat(y_ind,x_ind) = var_sum_mat(y_ind,x_ind)+firing_var_t(i);        
    end
end

var_map= var_sum_mat;
sqrt_rate_mat= sqrt(rate_mat);


%var_mean_mat= var_sum_mat./rate_mat;

% create rate mat
%var_map=var_map./time_mat;

 var_map=SmoothRateMat(var_map,parms);
 
 var_map= var_map./ sqrt_rate_mat;  
 
 
%var_map= var_map./rate_mat; %normalize by diving var by mean
% spike_mat_smooth =SmoothRateMat(spike_mat, parms);  % want to see spike_mat irrelevant to time 
 
 disp('');