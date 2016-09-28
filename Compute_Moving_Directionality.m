function [MD_rate_hist,MD_rayleigh_score,MD_rayleigh_angle,...
        max_hd_rate,pos_MD,spk_MD]=Compute_Moving_Directionality...
    (pos_t,pos_x,pos_y,spk_t,pos_vx,pos_vy,parms)
%this function recieves positional intput and a speed threshold
%outputs:
%rate_hist:histogram of MD rates normalized to time.
%rayscore:Rayleigh Score
%rayangle:Rayleigh angle

dt=median(diff(pos_t));

% the moving direction of the animal through out the trail
% pos_MD= atan2(diff(pos_y), diff(pos_x));
% 
% if isrow(pos_MD)
%    pos_MD=pos_MD';
% end
% pos_MD=[pos_MD; pos_MD(end)];
pos_MD= atan2(pos_vy,pos_vx);
spk_MD=interp1(pos_t,unwrap(pos_MD'),spk_t);
spk_MD=wrapToPi(spk_MD);

%% throw times when the rat was standing still
pos_speed=sqrt(pos_vx.^2+pos_vy.^2);
spk_speed=interp1(pos_t,pos_speed,spk_t);
threshold=2;
pos_MD=pos_MD(pos_speed>threshold);
spk_MD=spk_MD(spk_speed>threshold);
%%

n_bins=parms.num_of_direction_bins;
%ang_ax =0:2*pi/(n_bins-1):2*pi;
ang_ax =-pi:2*pi/(n_bins-1):pi;

hist_count_spk = hist(spk_MD,ang_ax);
hist_time = hist(pos_MD,ang_ax)*dt;

tmp_MD_rate_hist=hist_count_spk./hist_time;

%%convolution with hamming window
Win=hamming(10);
Win=Win/sum(Win);
tmp_MD_rate_hist=cconv(tmp_MD_rate_hist,Win');

%% SHAHAF PATCH
%this line is meant to solve a problem where the in line 22 the smoothing
%causes slightly negative values to apear.
tmp_MD_rate_hist(tmp_MD_rate_hist<0)=0;


MD_rate_hist=tmp_MD_rate_hist((length(Win)/2):length(tmp_MD_rate_hist)-(length(Win)/2));
[max_hd_rate,max_bin]=max(MD_rate_hist);

%% find rayleigh score & angle
norm_val = nansum(MD_rate_hist);
x_vec = nansum(cos(ang_ax).*MD_rate_hist);
y_vec = nansum(sin(ang_ax).*MD_rate_hist);
vec_len = sqrt(x_vec.^2+y_vec.^2);

MD_rayleigh_score = vec_len/norm_val;
MD_rayleigh_angle = atan2(y_vec,x_vec);

 %figure;bar(rad2deg(ang_ax),MD_rate_hist);
%     title('number of spikes fired in each direction');

%figure;bar(ang_ax,MD_MD_rate_hist);

disp('')