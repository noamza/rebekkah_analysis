function [rate_hist,rayleigh_score,rayleigh_angle]=ComputeMovingDirectionality...
    (pos_t,pos_x,pos_y,spk_t,parms)

[rate_hist,rayleigh_score,rayleigh_angle] = ComputeMovingDirectionalityWithSpeedThreshHold...
    (pos_t,pos_x,pos_y,spk_t,0,parms);


