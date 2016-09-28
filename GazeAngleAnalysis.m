parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\smooted gaze rate results- only one session included [for checking which wall tend sto look at]';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
    [size_x, size_y] = size(S.gaze_rate);
    
    top_gaze_rate= S.gaze_rate;
    top_gaze_rate(2:size_y,1:size_x)=nan;
    bottom_gaze_rate= S.gaze_rate;
    bottom_gaze_rate(1:size_y-1,1:size_x)=nan;
    left_gaze_rate= S.gaze_rate;
    left_gaze_rate(1:size_y,2:size_x)=nan;
    right_gaze_rate= S.gaze_rate;
    right_gaze_rate(1:size_y,1:size_x-1)=nan;
    
    [~,max_gaze_pt_top]= max(top_gaze_rate(:));
     [~,max_gaze_pt_bottom]= max(bottom_gaze_rate(:));
      [~,max_gaze_pt_left]= max(left_gaze_rate(:));
       [~,max_gaze_pt_right]= max(right_gaze_rate(:));
    
    [R,C] = ind2sub(size(S.gaze_rate),max_gaze_pt_top);
    max_gaze_pt_top=[];
    max_gaze_pt_top(1)=R;
    max_gaze_pt_top(2)=C;
    
    [R,C] = ind2sub(size(S.gaze_rate),max_gaze_pt_bottom);
    max_gaze_pt_bottom=[];
    max_gaze_pt_bottom(1)=R;
    max_gaze_pt_bottom(2)=C;
    
    [R,C] = ind2sub(size(S.gaze_rate),max_gaze_pt_left);
    max_gaze_pt_left=[];
    max_gaze_pt_left(1)=R;
    max_gaze_pt_left(2)=C;
   
    [R,C] = ind2sub(size(S.gaze_rate),max_gaze_pt_right);
    max_gaze_pt_right=[];
    max_gaze_pt_right(1)=R;
    max_gaze_pt_right(2)=C;
       
%     mid_pt= [size_x/2, size_y/2]; 
    
%     gaze_angle_t_b = atan2(abs(det([max_gaze_pt_bottom-mid_pt;max_gaze_pt_top-mid_pt])),dot(max_gaze_pt_bottom-mid_pt,max_gaze_pt_top-mid_pt));
%     gaze_angle_t_b= rad2deg(gaze_angle_t_b);
%     
%     gaze_angle_l_r = atan2(abs(det([max_gaze_pt_left-mid_pt;max_gaze_pt_right-mid_pt])),dot(max_gaze_pt_left-mid_pt,max_gaze_pt_right-mid_pt));
%     gaze_angle_l_r= rad2deg(gaze_angle_l_r);
   
pt_y= max_gaze_pt_top(1);
pt_x= max_gaze_pt_top(2);
m_y= max_gaze_pt_bottom(1);
m_x= max_gaze_pt_bottom(2);
        
        gaze_angle_t_b = abs(atan2(pt_y-m_y, pt_x-m_x) * 180/pi);
        
pt_y= max_gaze_pt_left(1);
pt_x= max_gaze_pt_left(2);
m_y= max_gaze_pt_right(1);
m_x= max_gaze_pt_right(2);

gaze_angle_l_r= abs(atan2(pt_y-m_y, pt_x-m_x) * 180/pi);

    gaze_angles(h,1)= gaze_angle_t_b;
    gaze_angles(h,2)= gaze_angle_l_r;
        
end