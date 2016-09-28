parms.dir_load_data = 'C:\Users\Dori\Desktop\Rebekkah data\data from non HD PF at least 8 cells\only cells with number of PF over 7 non HD'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% peak_pts_mat = zeros(100); 

for i=1:length(file_names)
    
    file_name = file_names{i};
    load(file_name);
    
    
    [size_x, size_y] = size(S.rate_mat);

    norm_max_index = []; 
    
    %% for actual data
%     
    norm_max_index(1,1) = S.max_index(1,1) / size_x;
    norm_max_index(1,2) = S.max_index(1,2) / size_y;
%     
   
    
    %% shuffles results
    
%     [size_x, size_y] = size(S.peak_zone_mat);
% 
%     peak_rates_list = S.peak_rates(randperm(length(S.peak_rates)));
% 
%     zone_num = find(peak_rates_list==max(peak_rates_list));
%     zone_num = zone_num(1,1);
%     
%     [row, col] = find(S.peak_zone_mat == zone_num); 
% 
%     shuffle_index = [];
%     shuffle_index(1,1) = row;
%     shuffle_index(1,2) = col;
%    
%     norm_max_index(1,1) = shuffle_index(1,1) / size_x;
%     norm_max_index(1,2) = shuffle_index(1,2) / size_y;

%%turn max peaks into matrix

%     pt1 = norm_max_index(1,2) * 100;
%     pt2 = norm_max_index(1,1) * 100;
% 
%     peak_pts_mat(int16(pt1), int16(pt2))= 1;

%% plot with x and y axis
    
    plot (norm_max_index(1,2), norm_max_index(1,1), 'x', 'MarkerSize', 15, 'color', 'r', 'LineWidth', 2);
    hold on;

end

%%image as matrix 

% imagesc(peak_pts_mat)
% plot(find(peak_pts_mat==1), 'x', 'MarkerSize', 15)



% xvalues= 0:0.5:5;
% hist(fano_factor(:), xvalues);
% xlim ([0 5]); 

disp(''); 