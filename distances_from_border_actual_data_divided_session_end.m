function distances_from_border_actual_data_divided_session_end
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    
    % get all files in folder
    % create file_list
    % create loop
    
    parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\examine first half and last half\results';
    
    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    distances = [];
    
    for i=1:length(file_names)
        
        file_name = file_names{i};
        load(file_name);
        
        [size_x, size_y] = size(S.rate_mat_end);
                
        max_peak_distance_right_left= [];
        max_peak_distance_top_bottom= [];
        max_peak_distance = []; 
 
      %% finds distances  
      
        
        h = 1;
        
        for cen = [1, size_x]
            for cen2 = 1:size_y
            max_peak_distance_top_bottom (1,h) = Distance(S.max_index_end(1, 1),S.max_index_end(1, 2), cen , cen2);
        
            h = h+1;
            end
        end
        
        h = 1;
        
        for cen = [1, size_y]
            for cen2= 1:size_x
            max_peak_distance_right_left(1,h)= Distance(S.max_index_end(1, 1),S.max_index_end(1, 2), cen2 , cen);
        
            h= h+1;
            end
        end
        
        max_peak_distance = union(max_peak_distance_top_bottom, max_peak_distance_right_left);
        
        max_peak_distance = min(max_peak_distance);
        
        [size_x, size_y] = size(S.rate_mat_end);
        %
        distances(i) = max_peak_distance/size_x;
        
        
    end
        
figure;
hist (distances(:));
 
save('distances PF 8 bin6 nonsmooth end', 'distances');

disp('');


