parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results 202 cells with adjusted zone mats [use to find PF sizes] [10 percent removed]';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
    location_inds = find(S.location == 1);
    
    %only look at examples where the max peak is touching one wall
    %(corner and central max peaks removed)
    
%    S.smallest_degree_wall = S.smallest_degree_Wall;
    
    if length(location_inds)==1
        max_peak_location = location_inds;
        
        
        if max_peak_location == 1 % if max peak located top wall
            if S.smallest_degree_wall == 1 | S.smallest_degree_wall== 2
                peak_angle_diff = 0; % location of max peak and min angle same
            elseif S.smallest_degree_wall == 3 | S.smallest_degree_wall== 4
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
            
        elseif max_peak_location == 2 % if max peak located bottom wall
            if S.smallest_degree_wall== 1 | S.smallest_degree_wall == 2
                peak_angle_diff = 0; % location of max peak and min angle same
           elseif S.smallest_degree_wall == 3 | S.smallest_degree_wall== 4
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
            
        elseif max_peak_location == 3 % if max peak located left wall
            if S.smallest_degree_wall == 3 | S.smallest_degree_wall == 4
                peak_angle_diff = 0; % location of max peak and min angle same
            elseif S.smallest_degree_wall == 1 | S.smallest_degree_wall== 2
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
            
        elseif max_peak_location == 4 % if max peak located right wall
            if S.smallest_degree_wall == 3 | S.smallest_degree_wall == 4
                peak_angle_diff = 0; % location of max peak and min angle same
             elseif S.smallest_degree_wall == 1 | S.smallest_degree_wall== 2
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
        end
        
    else
        
        peak_angle_diff = NaN;  % if max peak is in corner or center, ignore in analysis
    end
    
    same_or_diff(h) = peak_angle_diff;
    
end

figure; hist(same_or_diff(:));

save('smallest angle vs max peak location', 'same_or_diff');

disp('');