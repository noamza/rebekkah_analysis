parms.dir_load_data = 'C:\Users\Dori\Desktop\Rebekkah data\data from non HD PF at least 14 cells';
    
    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    for i=1:length(file_names)
        
        file_name = file_names{i};
        load(file_name);
        
        second_peak = [];
        second_peak = S.sorted_means(end-1);
       
        %shuffles results
        
        peak_rates_list = S.peak_rates(randperm(length(S.peak_rates)));
        
        zone_num = find(peak_rates_list==second_peak);
        
        % for actual data
        
%         zone_num = find(S.peak_rates==second_peak);
        
        %for both actual data and shuffled data
        
        [row, col] = find(S.peak_zone_mat == zone_num); 

        index = [];
        index(1,1) = row;
        index(1,2) = col;
   
        [size_x, size_y] = size(S.peak_zone_mat);
        
        norm_index(1,1) = index(1,1) / size_x;
        norm_index(1,2) = index(1,2) / size_y;
       
        plot (norm_index(1,2), norm_index(1,1), 'x', 'MarkerSize', 15, 'color', 'r', 'LineWidth', 2);
        hold on;
        
    end
        
disp('');         