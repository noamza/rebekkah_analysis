parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results with smoothed gaze mat and disregard when near border [0.05';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

    gaze_walls_0=[];
    gaze_walls_1=[];
    gaze_walls_2=[];
    gaze_walls_3=[];
    gaze_walls_4=[];
    

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
    gaze_wall=[];
    gaze_orient_diff=[];
    
    % find if location of maximum gazing correlates to orientation wall
    
    [size_x, size_y] = size(S.gaze_rate);
    [value, gaze_ind] = max(S.gaze_rate(:));
    
    [R,C] = ind2sub(size(S.gaze_rate),gaze_ind);
    gaze_ind(1)=R;
    gaze_ind(2)=C;
    
    if gaze_ind(1) >= 1 && gaze_ind(1) <=4
        gaze_wall = 1; %top wall
    elseif gaze_ind(1) <= size_x && gaze_ind(1) >= size_x-4
        gaze_wall = 2; %bottom wall
    elseif gaze_ind(2)>= 1 && gaze_ind(2)<= 4
        gaze_wall = 3; %left wall
    elseif gaze_ind(2)<= size_y && gaze_ind(2)>= size_y-4
        gaze_wall =4; %right wall
    else
        disp('ERROR ERROR')
    end
    
    
    
    if gaze_wall == 1 % if max peak located top wall
        if S.main_axis == 1
            gaze_orient_diff = 0; % location of gaze peak and orient axis same
        elseif S.main_axis==0
            gaze_orient_diff= 1; % location of max peak and min angle diff
        end
        
    elseif gaze_wall == 2 % if max peak located bottom wall
        if S.main_axis == 1
            gaze_orient_diff = 0; % location of max peak and min angle same
        elseif S.main_axis == 0
            gaze_orient_diff = 1; % location of max peak and min angle diff
        end
        
    elseif gaze_wall == 3 % if max peak located left wall
        if S.main_axis == 0
            gaze_orient_diff= 0; % location of max peak and min angle same
        elseif S.main_axis == 1
            gaze_orient_diff = 1; % location of max peak and min angle diff
        end
        
    elseif gaze_wall == 4 % if max peak located right wall
        if S.main_axis == 0
            gaze_orient_diff = 0; % location of max peak and min angle same
        elseif S.main_axis == 1
            gaze_orient_diff = 1; % location of max peak and min angle diff
        end
    end
    
    
    gaze_same_or_diff(h) = gaze_orient_diff;
    
  
    %gaze peak correlated to max peak? wall inds
    
      location_inds = find(S.location == 1);
    
    if length(location_inds)==1
        
     if location_inds == 4 % if max peak located top wall
        if gaze_wall == 1
            gaze_peak_diff = 0; % location of gaze peak and orient axis same
        else
            gaze_peak_diff= 1; % location of max peak and min angle diff
        end
        
    elseif location_inds==3% if max peak located bottom wall
        if gaze_wall == 2
            gaze_peak_diff = 0; % location of max peak and min angle same
        else
           gaze_peak_diff = 1; % location of max peak and min angle diff
        end
        
    elseif location_inds == 1 % if max peak located left wall
        if gaze_wall == 3
            gaze_peak_diff= 0; % location of max peak and min angle same
        else
            gaze_peak_diff = 1; % location of max peak and min angle diff
        end
        
    elseif location_inds == 2 % if max peak located right wall
        if gaze_wall == 4
            gaze_peak_diff = 0; % location of max peak and min angle same
        else
            gaze_peak_diff = 1; % location of max peak and min angle diff
        end
    end
        
    else
        
        gaze_peak_diff = NaN;  % if max peak is in corner or center, ignore in analysis
    end
    
    gaze_peak_diffs(h)=gaze_peak_diff;
    
    if S.axis_and_orient == 1
        gaze_walls_1(length(gaze_walls_1)+1)= gaze_wall;
    elseif S.axis_and_orient == 2
        gaze_walls_2(length(gaze_walls_2)+1) = gaze_wall;
    elseif S.axis_and_orient == 3
        gaze_walls_3(length(gaze_walls_3)+1)= gaze_wall;
    elseif S.axis_and_orient == 4
        gaze_walls_4(length(gaze_walls_4)+1)= gaze_wall;
    end
    
    gaze_walls(h)= gaze_wall;
    
end

figure; hist(gaze_same_or_diff(:)); %gaze and angle orientation [0.5 chance probability]
figure; hist(gaze_peak_diffs(:));   %gaze and wall of peak [0.25 chance probabilty]
figure; subplot(1,5,1); hist(gaze_walls_1(:));
subplot(1,5,2); hist(gaze_walls_2(:));
subplot(1,5,3); hist(gaze_walls_3(:));
subplot(1,5,4); hist(gaze_walls_4(:));
subplot(1,5,5); hist(gaze_walls);

save('gaze vs orientation same or diff', 'gaze_same_or_diff', 'gaze_peak_diffs', 'gaze_walls');

disp('');