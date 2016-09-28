function second_max_peaks_near_borders 
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


%for k= 1:500
    
    % get all files in folder
    % create file_list
    % create loop
    
    parms.dir_load_data = 'N:\users\rebekkah\bin size 6 nonsmooth\results updated';
    
    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    distances = [];
    
    for i=1:length(file_names)
        
        file_name = file_names{i};
        load(file_name);
                
        second_peak = S.sorted_means(end-1);
       
        % shuffles results
        
%         peak_rates_list = S.peak_rates(randperm(length(S.peak_rates)));
%         
%         zone_num = find(peak_rates_list==second_peak);
        
        % for actual data
        
        zone_num = find(S.peak_rates==second_peak);
        
        
      %  for both actual data and shuffled data   
          
       zone_num = zone_num(1); % if same firing rate at two locations, just take first for analysis
      
        [row, col] = find(S.peak_zone_mat == zone_num); 

        index = [];
        index(1,1) = row;
        index(1,2) = col;
   
        [size_x, size_y] = size(S.peak_zone_mat);
            
        second_max_peak_distance_right_left= [];
        second_max_peak_distance_top_bottom= [];
        second_max_peak_distance= [];
        
        %% finds distance for shuffled data
        
%         h = 1;
% %         
%         for cen = [1, size_x]
%             for cen2 = 1:size_y
%                 second_max_peak_distance_top_bottom (1,h) = Distance(index(1, 1),index(1, 2), cen , cen2);
%                 
%                 h = h+1;
%             end
%         end
%         
%         h = 1;
%         
%         for cen = [1, size_y]
%             for cen2= 1:size_x
%                 second_max_peak_distance_right_left(1,h)= Distance(index(1, 1),index(1, 2), cen2 , cen);
%                 
%                 h= h+1;
%             end
%         end
%         
%         second_max_peak_distance = union(second_max_peak_distance_top_bottom, second_max_peak_distance_right_left);
%         
%         second_max_peak_distance = min(second_max_peak_distance);
        
        
        %normalizing of distances of border to peak firing field
        
%         [size_x, size_y] = size(S.rate_mat);
%         
%         distances(i) = second_max_peak_distance/size_x;
        
        
        
%         %% distances for actual data
%         
        h = 1;
        
        for cen = [1, size_x]
            for cen2 = 1:size_y
            second_max_peak_distance_top_bottom (1,h) = Distance(index(1, 1),index(1, 2), cen , cen2);
        
            h = h+1;
            end
        end
        
        h = 1;
        
        for cen = [1, size_y]
            for cen2= 1:size_x
            second_max_peak_distance_right_left(1,h)= Distance(index(1, 1),index(1, 2), cen2 , cen);
        
            h= h+1;
            end
        end
        
        second_max_peak_distance = union(second_max_peak_distance_top_bottom, second_max_peak_distance_right_left);
        
        second_max_peak_distance = min(second_max_peak_distance);
        
        [size_x, size_y] = size(S.rate_mat);
       
        distances(i) = second_max_peak_distance/size_x;
        
        
    end
    
    
   save('second max peak distances from border stats', 'distances');
    
    %
    
    
%     hist(distances(:));
    %
    % %counts number of cell with highest peaking PF near border
    %
%     percent_near_border = sum(distances < 0.05) ;
%     
%     percent2 = sum (distances < 0.1);
%     
%     percent3 = sum(distances < 0.15); 
%     
%     mean_of_set= mean(distances);
%     % percent_at_zero = sum(distances ==0);
%     %
%     border_occurances(k) = percent_near_border; % less than 0.05
%     
%     border2 (k) = percent2;     % less than 0.1
%     
%     border3 (k) = percent3;     %less than 0.15
%     
%     border_means (k)= mean_of_set;
%     % border_at_zero(k) = percent_at_zero;
%end

% figure; 
% n =1;
% m=3;
% 
% subplot(n,m,1)
% hist (border_occurances(:));
% 
% subplot(n,m,2)
% hist(border2(:));
% 
% subplot(n,m,3)
% hist(border3(:));
% 
% save('PF8 bin6 nonsmooth second peak shuffled AGAIN', 'border_occurances', 'border2', 'border3', 'border_means')
% 
% disp('');


