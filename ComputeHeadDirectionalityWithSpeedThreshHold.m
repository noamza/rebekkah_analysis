function [rate_hist,ray_score,ray_angle]=ComputeHeadDirectionalityWithSpeedThreshHold...
    (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_t,speed_thresh_hold,parms)
    %(pos_x,pos_y,pos_x2,pos_y2,spk_x,spk_y,spk_x2,spk_y2,parms)

dt = pos_t(2) - pos_t(1);
speed = sqrt(diff(pos_x).^2 + diff(pos_y).^2) / dt;
    
% the head direction of the animal through out the trail
time_phi = wrapTo2Pi(atan2(pos_y2-pos_y,pos_x2-pos_x));
spk_phi = interp1(pos_t,time_phi,spk_t);
spk_speed = interp1(pos_t,[0 speed'],spk_t);
% cut all the data under the thresh hold
time_phi = time_phi(find(speed>=speed_thresh_hold));
count_phi = spk_phi(find(spk_speed>=speed_thresh_hold));

step=deg2rad(3);
ang_ax = 0:step:2*pi;
count = hist(count_phi,ang_ax);
time = hist(time_phi,ang_ax)*dt;

[rate_hist,ray_score,ray_angle] = CalculateRayliehAndRateHist(count,time,ang_ax);

%figure;bar(rad2deg(ang_ax),rate_hist);title('number of spikes fired in each direction');

