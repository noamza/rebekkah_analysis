
parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA improved maybe';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
  % [size_x, size_y]= size(rate_mat);

    %envir_x= [1:size_x ones(1,size_y-1)*(size_x) size_x-1:-1:1 ones(1,size_y-2)];

    %envir_y= [ones(1,size_x) [2:size_y] ones(1,size_x-1)*size_y [size_y-1:-1:2]];

   figure; hist(S.hist_count_gaze); 
    
    
    % find most looked at location
%     binranges_corner_4= [-135 -45 45 135];
%     gaze_location_corner= histc(location_of_gaze_in_degrees,binranges_corner_4);
%     
%     binranges_wall_4= [-180 -90 0 90 180];
%     gaze_location_wall= histc(location_of_gaze_in_degrees,binranges_wall_4);
%     
%     binranges_8_corner= [-180 -126 -81 -36 9 67.5 112.5 157.5];
%     gaze_location_8_corner= histc(location_of_gaze_in_degrees, binranges_8_corner);
%     
%     binranges_8_wall= [-157.5 -112.5 -67.5 -22.5 22.5 67.5 112.5 157.5];
%     gaze_location_8_wall= histc(location_of_gaze_in_degrees, binranges_8_wall);
%     
% %     binranges_top_bottom= [0 180];
% %     gaze_location_t_b= histc(location_of_gaze_in_degrees, binranges_top_bottom);
% %     
% %     binranges_left_right= [-90 90];
% %     gaze_location_l_r= histc(location_of_gaze_in_degrees, binranges_left_right);
%     
%     gaze_direction_corner= binranges_4(gaze_location_corner==max(gaze_location_corner));
%     gaze_direction_rough= binranges_4(gaze_location_rough==max(gaze_location_rough));
%     gaze_direction_specific= binranges_8(gaze_location_specific==max(gaze_location_specific));

 %   max_index_degrees= atan2(S.max_index(2)-middle_pt_y, S.max_index(1)-middle_pt_x);
  %  max_index_degrees= radtodeg(max_index_degrees)
%max_index_degrees= max_index_degrees * 360/(2*pi) 
   
    

%     if max_index_degrees > gaze_direction_rough - 45 && max_index_degrees < gaze_direction_rough+45
%         gaze_and_pivot_pt_rough(h)= 0; %same
%     else
%         gaze_and_pivot_pt_rough(h)= 2; %diff
%     end
%     
%     if max_index_degrees > gaze_direction_specific - 22.5 && max_index_degrees < gaze_direction_specific+22.5
%         gaze_and_pivot_pt_specific(h)= 0; %same
%     else
%         gaze_and_pivot_pt_specific(h)= 2; %diff
%     end
    

% score based on how close the gaze is 

% if location_inds_len==2 %if max-peak in corner
%     if max_index_degrees > gaze_direction_specific - 22.5 && max_index_degrees < gaze_direction_specific+22.5
%         gaze_and_pivot_pt(h)= 0; %same to 45 degrees 
%     elseif max_index_degrees > gaze_direction_rough - 45 && max_index_degrees < gaze_direction_rough+45
%         gaze_and_pivot_pt(h)= 1; %same to 90 degrees  
%     else
%         gaze_and_pivot_pt(h)= 2; %diff
%     end
%     
%    


%     if max_index_degrees > gaze_direction_t_b - 90 && max_index_degrees < gaze_direction_t_b+90
%         gaze_and_pivot_pt_t_b(h)= 0; %same
%     else
%         gaze_and_pivot_pt_t_b(h)= 2; %diff
%     end
%     
%     if max_index_degrees > gaze_direction_l_r - 90 && max_index_degrees < gaze_direction_l_r+90
%         gaze_and_pivot_pt_l_r(h)= 0; %same
%     else
%         gaze_and_pivot_pt_l_r(h)= 2; %diff
%     end
    
    disp('');
    
end

%figure; hist(gaze_and_pivot_pt);
%figure; hist(gaze_and_pivot_pt_specific);