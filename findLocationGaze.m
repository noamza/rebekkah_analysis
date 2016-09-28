
function [location_of_gaze_x, location_of_gaze_y, hist_count_gaze, gaze_rate, gaze_spike_rate]=findLocationGaze(pos_HD,parms,rate_mat,pos_t, pos_mean_y, pos_mean_x, spk_t)

  location_of_gaze_x=[];
  location_of_gaze_y=[];
  hist_count_gaze=[];

pos_MD=wrapToPi(atan2(diff(pos_mean_y), diff(pos_mean_x)));
pos_MD=[pos_MD;pos_MD(end)];

 pos_HD=wrapTo180(rad2deg(pos_HD));
 pos_MD=wrapTo180(rad2deg(pos_MD));
 
 HD_MD=wrapTo180(pos_HD-pos_MD);
 
 HD_MD(HD_MD>-180)=HD_MD(HD_MD>-180)-180;
 HD_MD(HD_MD<180)=HD_MD(HD_MD<180)+180;

 n_bins=parms.num_of_direction_bins;
 %ang_ax =-pi:2*pi/(n_bins-1):pi;

[diff_HD_MD_hist,ang_ax] = hist(HD_MD,n_bins);
[max_val,max_ind] = max(diff_HD_MD_hist);
 HD_subtract_angle = ang_ax(max_ind);
  pos_HD_corrected = wrapTo180(pos_HD - HD_subtract_angle);
 
 pos_HD_corrected=deg2rad(pos_HD_corrected);
 HD_subtract_angle=deg2rad(HD_subtract_angle);
  
 pos_HD_corrected=pos_HD_corrected';
 
 %% disregard when the rat is near the border
  
%  for h= 1:length(pos_t)
%  if max(abs(pos_mean_x(h)),abs(pos_mean_y(h))) >= sqrt(0.5)*max(pos_mean_x) %removed max
%      plot(pos_mean_x(h), pos_mean_y(h))
%     % pos_mean_x(h)= nan;
%     % pos_mean_y(h)= nan; 
%  end
% end
%  



%% find location rat is looking towards

% convert position coordinates to rate mat coordinates

[size_x, size_y]= size(rate_mat);

bin_size=3;

pos_x_converted= pos_mean_x;
pos_x_converted=pos_x_converted';
nan_inds=find(isnan(pos_mean_x));
pos_x_converted=pos_mean_x/bin_size + size_x/2; %remove round
pos_x_converted(nan_inds)= nan;  
pos_x_converted=pos_x_converted';

pos_y_converted= pos_mean_y;
pos_y_converted=pos_y_converted';
nan_inds=find(isnan(pos_mean_y));
pos_y_converted= pos_mean_y/bin_size + size_y/2; %remove round
pos_y_converted(nan_inds)= nan; 
pos_y_converted=pos_y_converted';


%%% disregard 10% of border positions 
for h= 1:length(pos_t)
    if pos_x_converted(h) >= 0 & pos_x_converted(h) <= 0.15*size_y ...
        | pos_x_converted(h) >= 0.85*size_y & pos_x_converted(h) <= size_y ...
        | pos_y_converted(h) >= 0 & pos_y_converted(h) <= 0.15*size_x...
        | pos_y_converted(h) >= 0.85*size_x & pos_y_converted(h)<= size_x
        pos_x_converted(h) = nan;
        pos_y_converted(h) = nan;
    end
end

envir_x= [1:size_x ones(1,size_y-1)*(size_x) size_x-1:-1:1 ones(1,size_y-2)];

envir_y= [ones(1,size_x) [2:size_y] ones(1,size_x-1)*size_y [size_y-1:-1:2]];

for h=1:length(pos_t);

    direction_in_environment= atan2(pos_y_converted(h)-envir_y, pos_x_converted(h)-envir_x);
     
    if ~isnan(pos_x_converted(h)) && ~isnan(pos_y_converted(h))
        
        diff_vec = direction_in_environment - pos_HD_corrected(h);
        diff_vec(diff_vec > pi) = diff_vec(diff_vec > pi) - pi;
        diff_vec(diff_vec < -pi) = diff_vec(diff_vec < -pi) + pi;
        [value, ind]= nanmin(abs(diff_vec));
        
        location_of_gaze_x(h)= envir_x(ind);    %gaze location at each time stamp (size of pos_t)
        location_of_gaze_y(h)= envir_y(ind);
        hist_count_gaze(h)= ind;                %index of gaze location at each time stamp (size of pos_t)
        
    else
        location_of_gaze_x(h)= NaN;
        location_of_gaze_y(h)= NaN;
        hist_count_gaze(h)= NaN;
    end
end



 %HD_MD_corrected_new =wrapTo180(pos_HD_corrected-pos_MD);
%% 90 offset correction, not currently used: 
%  HD_MD_corrected =wrapTo180(pos_HD_corrected-pos_MD);
%  exchange_inds = HD_MD_corrected > 90 | HD_MD_corrected < -90;
%  pos_HD_corrected(exchange_inds) = wrapTo180(pos_HD_corrected(exchange_inds)+180);
%%
gaze_rate = nan*zeros(size(rate_mat));
for ind = 1:length(envir_x)
   num(ind) = length(find(hist_count_gaze == ind));
   gaze_rate(envir_x(ind),envir_y(ind)) = num(ind);
end
% corner_inds = find(hist_count_gaze == 1);
% plot(pos_mean_y(corner_inds),pos_mean_x(corner_inds),'.');
% figure;
% imagesc(gaze_rate)

%% finds spikes rate map of gaze location


gaze_location_at_each_spike = interp1(pos_t,hist_count_gaze,spk_t,'nearest');  %gives the gaze location at each spike time (size of spk_t)
        
gaze_spike_rate = nan*zeros(size(rate_mat));
for ind = 1:length(envir_x)
   num(ind) = length(find(gaze_location_at_each_spike == ind));     %finds number of times at each border position that spikes 
   gaze_spike_rate(envir_x(ind),envir_y(ind)) = num(ind);
end


disp('')