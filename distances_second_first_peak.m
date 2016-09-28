parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results sorted PF 93 cells'; 

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=6;
%parms.sigma = 3;
  
dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count=0;

figure;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};    
    load(file_name); 
    
    %to shuffle add this
    
     second_peak= S.sorted_means(end-1);
    
     peak_rates_list = S.peak_rates(randperm(length(S.peak_rates)));    
     zone_num = find(peak_rates_list==second_peak);
    
     [row, col] = find(S.peak_zone_mat == zone_num);
        
        second_index(1) = row;
        second_index(2) = col;
     
        peak= S.sorted_means(end);   
     zone_num = find(peak_rates_list==peak);
    
     [row, col] = find(S.peak_zone_mat == zone_num);
        
        max_index(1) = row;
        max_index(2) = col;
        
%     if S.sorted_means(end-1) >= S.sorted_means(end)*0.95
%         
%         % find distance between the two points
%         
%         % find index of second point
%         %find value in zone_mat
%         %use index to find number_zone
%         %use number zone to find point
%         
%         second_peak= S.sorted_means(end-1);
%         zone_num = find(S.peak_rates==second_peak);
%         [row, col] = find(S.peak_zone_mat == zone_num);
%         
%         second_index(1) = row;
%         second_index(2) = col;
%    
%         [size_x, size_y] = size(S.peak_zone_mat);
%         
%         norm_2_index(1) = second_index(1) / size_x;
%         norm_2_index(2) = second_index(2) / size_y;
%         
%         distance_of_points= Distance(S.norm_max_index(1), norm_2_index(1), S.norm_max_index(2), norm_2_index(2));
%         
%         count=count+1;
%         
%           distances_of_points(count)= distance_of_points;       
%     end
%   
    
%for actual data

%      second_peak= S.sorted_means(end-1);
%         zone_num = find(S.peak_rates==second_peak);
%         [row, col] = find(S.peak_zone_mat == zone_num);
%         
%         second_index(1) = row;
%         second_index(2) = col;
   
%for both

        [size_x, size_y] = size(S.peak_zone_mat);
        
        norm_2_index(1) = second_index(1) / size_x;
        norm_2_index(2) = second_index(2) / size_y;
        
        max_index(1)= max_index(1)/size_x;
        max_index(2)= max_index(2)/size_y;
        
        S.norm_max_index=max_index;
        
    second= (S.sorted_means(end-1) - mean(S.sorted_means))/std(S.sorted_means);
    first= (S.sorted_means(end) - mean(S.sorted_means))/std(S.sorted_means);
    
    PF_radius= S.PF_radius/size_x;
    
    ratio_second_first(i)= S.sorted_means(end-1)/S.sorted_means(end);
    distances(i)= Distance(S.norm_max_index(1), norm_2_index(1), S.norm_max_index(2), norm_2_index(2));
    distances(i)= distances(i)/PF_radius; %divide distance between points by grid spacing
    
    
    cell_nums(i)= S.i;
    
    
    if ratio_second_first(i) >= 0.9
        count=count+1;
        high_ratio_second_first(count)= ratio_second_first(i);
        high_ratio_PF_distances(count)= distances(i);
        high_ratio_cell_nums(count)= cell_nums(i);
    end
    
%     if S.number_of_PF >=10
%         count=count+1;
%         manyPF_ratio_second_first(count)= S.sorted_means(end-1)/S.sorted_means(end);
%         manyPF_distances(count)= Distance(S.norm_max_index(1), norm_2_index(1), S.norm_max_index(2), norm_2_index(2));
%     end
%   
% notes: 0.76- not considered high second to first max peak
%      0.84, 0.86- kind of high, not close together in the example
end

figure;
scatter(ratio_second_first, distances, 'o'); hold on;

%figure;
%scatter(manyPF_distances,manyPF_ratio_second_first, 'o'); hold on;

save('ratios and distances', 'ratio_second_first', 'distances', 'cell_nums');
save('high ratios', 'high_ratio_second_first', 'high_ratio_PF_distances', 'high_ratio_cell_nums');
   