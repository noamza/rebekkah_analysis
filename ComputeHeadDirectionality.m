function [rate_hist,rayleigh_score,rayleigh_angle]=ComputeHeadDirectionality...
    (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_t,parms)

    [rate_hist,rayleigh_score,rayleigh_angle] = ComputeHeadDirectionalityWithSpeedThreshHold...
        (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_t,0,parms);


