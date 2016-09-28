function [return_hist]=ComputePrecentOfMovingDirectionalityInHeadDirectionality...
    (pos_x,pos_y,pos_x2,pos_y2)
     % Assumption: return values between -pi...pi

% calculate the location of the rat by both leds
pos_average_x = (pos_x + pos_x2)/2;
pos_average_y = (pos_y + pos_y2)/2;

%% direction of rat movment at each moment
time_MD = wrapTo2Pi(atan2(diff(pos_average_y), diff(pos_average_x)));
if (size(time_MD,1) == 1)
   time_MD = time_MD'; 
end
time_MD_fix = [NaN; time_MD];

%% the direction of the rat looking at each moment
time_HD = wrapTo2Pi(atan2(pos_y2-pos_y,pos_x2-pos_x));
if (size(time_HD,1) == 1)
   time_HD = time_HD'; 
end
return_hist=time_MD_fix - time_HD;

% return values between -pi...pi
% notice: this is the same as using wrapToPi()
return_hist(return_hist>pi)=return_hist(return_hist>pi) - 2*pi;
return_hist(return_hist<-pi)=return_hist(return_hist<-pi) + 2*pi;

