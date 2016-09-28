function [pos_vx,pos_vy]=Compute_velocity_XY(mean_pos_x,mean_pos_y,pos_t,...
        parms)
   
    %convert to kate jefferys coordinates (about x5 cm) -rebekkah
    mean_pos_x=(mean_pos_x/5);
    mean_pos_y=(mean_pos_y/5);
       
   smoothing_parameter_for_velocity=parms.smoothing_factor_for_velocity;
    
    mean_pos_x_smooth = csaps(1:length(pos_t),...
             mean_pos_x, smoothing_parameter_for_velocity);
    mean_pos_y_smooth = csaps( 1:length(pos_t),...
            mean_pos_y, smoothing_parameter_for_velocity);

mean_pos_x_smooth=mean_pos_x;
mean_pos_y_smooth=mean_pos_y;
      
        % Compute the velocity via a derivative of the pp-form of the spline (fnder.m):
         tmp_pos_vx = fnval( fnder(mean_pos_x_smooth ),...
             1:length(pos_t) ) ;
         tmp_pos_vy = fnval( fnder(mean_pos_y_smooth),...
             1:length(pos_t) ) ;
             
        % Translate Velocity from pixels/frame to pixels/sec:
        delta_t_frame = mean( diff(pos_t));
        pos_vx =tmp_pos_vx / delta_t_frame; %* 10^6 ;
        pos_vy =tmp_pos_vy / delta_t_frame; %* 10^6 ;
        
    
disp('')