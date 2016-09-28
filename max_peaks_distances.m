function max_peaks_near_borders 
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


 for k= 1:500
    
    % get all files in folder
    % create file_list
    % create loop
    
    parms.dir_load_data = 'C:\Users\Dori\Desktop\bin size 6 nonsmooth\results';
    
    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    distances = [];
    
    for i=1:length(file_names)
        
        file_name = file_names{i};
        A = load(file_name);
               
        %%% shuffles results
        
        [size_x, size_y] = size(A.S.peak_zone_mat);
        
        peak_rates_list=[];
        
        peak_rates_list = A.S.peak_rates(randperm(length(A.S.peak_rates)));
        
        zone_num = find(peak_rates_list==max(peak_rates_list));
        zone_num = zone_num(1,1);
        
        
        [row, col] = find(A.S.peak_zone_mat == zone_num);
        
        shuffle_index = [];
        
        shuffle_index(1,1) = row;
        shuffle_index(1,2) = col;
        
        % %%% plots normalized location of highest peak firing rate for each cell
        %
        [size_x, size_y] = size(A.S.rate_mat);
        
        norm_max_index = [];
        
        %%% for shuffled data
        % %
        norm_max_index(1,1) = shuffle_index(1,1) / size_x;
        
        norm_max_index(1,2) = shuffle_index(1,2) / size_y;
        
        
        %% finds distance for shuffled data
        
        max_peak_distance_top_bottom=[];
        max_peak_distance_right_left=[];
        max_peak_distance=[];
        
        h = 1;
        
        for cen = [1, size_x]
            for cen2 = 1:size_y
                max_peak_distance_top_bottom (1,h) = Distance(shuffle_index(1, 1),shuffle_index(1, 2), cen , cen2);
                
                h = h+1;
            end
        end
        
        h = 1;
        
        for cen = [1, size_y]
            for cen2= 1:size_x
                max_peak_distance_right_left(1,h)= Distance(shuffle_index(1, 1),shuffle_index(1, 2), cen2 , cen);
                
                h= h+1;
            end
        end
        
        max_peak_distance = union(max_peak_distance_top_bottom, max_peak_distance_right_left);
        
        max_peak_distance = min(max_peak_distance);
        
        
       % normalizing of distances of border to peak firing field
        
        [size_x, size_y] = size(A.S.rate_mat);
        
        distances(i) = max_peak_distance/size_x;
        
   

    end
    
    
    %  hist(distances(:));
    %
    % %counts number of cell with highest peaking PF near border
    %
    percent_at_zero = sum(distances == 0);
    
    percent_near_border = sum(distances < 0.05) ;
    
    percent2 = sum (distances < 0.1);
    
    percent3 = sum(distances < 0.15); 
    %
    border_at_zero(k) = percent_at_zero;
    
    border_occurances(k) = percent_near_border;
    
    border2 (k) = percent2;
    
    border3 (k) = percent3;
 end

figure;
n =1;
m=3;

subplot(n,m,1)
hist (border_occurances(:));

subplot(n,m,2)
hist(border2(:)); hold on;
x= ones(1,140)*56;
y=0:139;
plot(x,y, 'color', 'r', 'linewidth', 2);

subplot(n,m,3)
hist(border3(:)); hold on;
x= ones(1,140)*67;
y=0:139;
plot(x,y, 'color', 'r', 'linewidth', 2);

save('PF8 bin6 nonsmooth shuffled', 'border_at_zero', 'border_occurances', 'border2', 'border3')

disp('');


