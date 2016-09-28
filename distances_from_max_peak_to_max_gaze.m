function distances_from_border_actual_data
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    
    % get all files in folder
    % create file_list
    % create loop
    
    parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results with smoothed gaze mat and disregard when near border [0.05';
    
    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    for i=1:length(file_names)
        file_name = file_names{i};
        load(file_name);
        
        [size_x, size_y] = size(S.gaze_spike_rate);
        [value, gaze_spike_ind] = max(S.gaze_spike_rate(:));
             
        [R,C] = ind2sub(size(S.gaze_spike_rate),gaze_spike_ind);
        
        norm_max_spike_gaze(1)= R/size_x;
        norm_max_spike_gaze(2)= C/size_y;
        
 
        [size_x, size_y] = size(S.gaze_rate);
        [value, gaze_value_ind] = max(S.gaze_rate(:));
              
        [R1,C1] = ind2sub(size(S.gaze_rate),gaze_value_ind);
        
        norm_max_gaze(1)= R1/size_x;
        norm_max_gaze(2)= C1/size_y;
        
        
      %% finds distances of max peak to max gaze 
      
      max_peak_max_gaze_distance = Distance(S.norm_max_index(1, 1),S.norm_max_index(1, 2), norm_max_gaze(1) , norm_max_gaze(2)); 
      
       max_peak_max_spike_gaze_distance = Distance(S.norm_max_index(1, 1),S.norm_max_index(1, 2), norm_max_spike_gaze(1) , norm_max_spike_gaze(2));   
      
      all_distances_peak_gaze(i) = max_peak_max_gaze_distance;
      
       all_distances_peak_spike_gaze(i) = max_peak_max_spike_gaze_distance;
      
      
        
    
        
        
    end
        
figure;
hist ( all_distances_peak_gaze(:));
figure;
hist ( all_distances_peak_spike_gaze(:)); 

save('distances of max peaks to max gaze peaks', 'all_distances_peak_gaze', 'all_distances_peak_spike_gaze');

disp('');


