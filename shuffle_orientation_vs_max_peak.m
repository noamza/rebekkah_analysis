
parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data with max peaks touching walls';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for k=1:1000

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);

[size_x,size_y]= size(S.zone_mat);

number_zone_mat= S.number_zone_mat;

%initialize everything

zones_that_touch_wall=[];
zones_that_touch_wall2=[];
zones_that_touch_wall3=[];
zones_that_touch_wall4=[];

zones_that_touch_wall_top_bottom=[];
zones_that_touch_wall_right_left=[];

corner_zones=[];
indices=[];

%finds zone numbers that touch walls 

[x,y, zones_that_touch_wall] = find(number_zone_mat(1,:));
[x,y, zones_that_touch_wall2]= find(number_zone_mat(:,size_y));
[x,y, zones_that_touch_wall3]= find(number_zone_mat(:,1));
[x,y, zones_that_touch_wall4]= find(number_zone_mat(size_x, :));

zones_that_touch_wall_top_bottom = union(zones_that_touch_wall,zones_that_touch_wall4);
zones_that_touch_wall_right_left= union(zones_that_touch_wall2,zones_that_touch_wall3);
zones_that_touch_wall_right_left = zones_that_touch_wall_right_left';

corner_zones= intersect(zones_that_touch_wall_top_bottom,zones_that_touch_wall_right_left);

zones_that_touch_wall = union(zones_that_touch_wall_top_bottom, zones_that_touch_wall_right_left);
    
    if ~isempty(corner_zones)
        if length(corner_zones) == 2
        zones_that_touch_wall(find(zones_that_touch_wall==corner_zones(1)))= []; %remove corner zones
        zones_that_touch_wall(find(zones_that_touch_wall==corner_zones(2)))= [];
        else
        zones_that_touch_wall(find(zones_that_touch_wall==corner_zones))= []; 
        end
    end
    
zones_wall_len= length(zones_that_touch_wall);

% randomly chooses one of the wall-touching-zones to be max 

rand_num = randi(zones_wall_len);

new_max_zone = zones_that_touch_wall(rand_num);

%finds location of this zone
[size_x, size_y]= size(S.zone_mat);

number_zone_mat=S.number_zone_mat;

top = unique(number_zone_mat(1,:));
top(find(top==0))= [];

left = unique(number_zone_mat(:,1));
left(find(left==0))= [];

right = unique(number_zone_mat(:,size_y));
right(find(right==0))= [];

bottom = unique(number_zone_mat(size_x,:)) ;  
bottom(find(bottom==0))= [];

location(1) = any(top==new_max_zone);
location(2) = any(bottom==new_max_zone);
location(3) = any(left==new_max_zone);
location(4) = any(right==new_max_zone);

location_inds = find(location == 1);

%% taken from individual_orientation_analysis

if length(location_inds)==1
        max_peak_location = location_inds;
        
        
        if max_peak_location == 1 % if max peak located top wall
            if S.smallest_degree_wall == 1 | S.smallest_degree_wall == 3
                peak_angle_diff = 0; % location of max peak and min angle same
            else
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
            
        elseif max_peak_location == 2 % if max peak located bottom wall
            if S.smallest_degree_wall == 2 | S.smallest_degree_wall == 4
                peak_angle_diff = 0; % location of max peak and min angle same
            else
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
            
        elseif max_peak_location == 3 % if max peak located left wall
            if S.smallest_degree_wall == 1 | S.smallest_degree_wall == 3
                peak_angle_diff = 0; % location of max peak and min angle same
            else
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
            
        elseif max_peak_location == 4 % if max peak located right wall
            if S.smallest_degree_wall == 2 | S.smallest_degree_wall == 4
                peak_angle_diff = 0; % location of max peak and min angle same
            else
                peak_angle_diff = 1; % location of max peak and min angle diff
            end
        end
        
    else
        
        peak_angle_diff = NaN;  % if max peak is in corner or center, ignore in analysis
    end
    
    same_or_diff(h) = peak_angle_diff;
    
    disp('');
    
end  

number_of_opposite(k)=sum(same_or_diff==1);

end

save('orientation vs max peak shuffle 1000x', 'number_of_opposite');

figure; hist(number_of_opposite(:));
hold on;
x= ones(1,140)*28;
y=0:139;
plot(x,y, 'color', 'r', 'linewidth', 2);

disp('');