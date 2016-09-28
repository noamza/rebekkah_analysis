parms.dir_load_data1 = 'N:\users\rebekkah\final data smoothed\results with fitted box arena';
parms.dir_load_data2 = 'N:\users\rebekkah\final data smoothed\major orientation results';

dir_name= parms.dir_load_data1;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count =1;
count2=1;

for h =1:length(file_names)
    cd(parms.dir_load_data1);
    file_name = file_names{h};
    dat1 = load(file_name);
    Cell1= dat1.S;
    
    cd(parms.dir_load_data2);
    dat2 = load(file_name);
    Cell2= dat2.S;
    
    max_pts_inside= find(Cell1.phi_of_arena.percent_inside_final(:,1)==max(Cell1.phi_of_arena.percent_inside_final(:,1)));
    
    arena_angle_changes= Cell1.phi_of_arena.percent_inside_final(max_pts_inside,5);
    
    angle_changes(h)= arena_angle_changes(1);
    
    if Cell2.wall==1 | Cell2.wall==4
        new_angle= Cell2.smallest_degree + arena_angle_changes(1);
    elseif Cell2.wall==2 | Cell2.wall==3
        new_angle= Cell2.smallest_degree - arena_angle_changes(1);
    end
    
%     location_inds = find(Cell2.location == 1);
%     
%     if new_angle >6.4 & new_angle <8.4
%         if length(location_inds)==1
%             max_peak_location = location_inds;
%             
%             if max_peak_location == 1 % if max peak located top wall
%                 if Cell2.smallest_degree_wall == 1 | Cell2.smallest_degree_wall == 3
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%                 
%             elseif max_peak_location == 2 % if max peak located bottom wall
%                 if Cell2.smallest_degree_wall == 2 | Cell2.smallest_degree_wall == 4
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%                 
%             elseif max_peak_location == 3 % if max peak located left wall
%                 if Cell2.smallest_degree_wall == 1 | Cell2.smallest_degree_wall == 3
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%                 
%             elseif max_peak_location == 4 % if max peak located right wall
%                 if Cell2.smallest_degree_wall == 2 | Cell2.smallest_degree_wall == 4
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%             end
%             
%         else
%             
%             peak_angle_diff = NaN;  % if max peak is in corner or center, ignore in analysis
%         end
%         
%         same_or_not_7_degree(count)= peak_angle_diff;
%         
%         count=count+1;
%     end
%     
%     if new_angle >6 & new_angle <8.8
%         if length(location_inds)==1
%             max_peak_location = location_inds;
%             
%             
%             if max_peak_location == 1 % if max peak located top wall
%                 if Cell2.smallest_degree_wall == 1 | Cell2.smallest_degree_wall == 3
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%                 
%             elseif max_peak_location == 2 % if max peak located bottom wall
%                 if Cell2.smallest_degree_wall == 2 | Cell2.smallest_degree_wall == 4
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%                 
%             elseif max_peak_location == 3 % if max peak located left wall
%                 if Cell2.smallest_degree_wall == 1 | Cell2.smallest_degree_wall == 3
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%                 
%             elseif max_peak_location == 4 % if max peak located right wall
%                 if Cell2.smallest_degree_wall == 2 | Cell2.smallest_degree_wall == 4
%                     peak_angle_diff = 0; % location of max peak and min angle same
%                 else
%                     peak_angle_diff = 1; % location of max peak and min angle diff
%                 end
%             end
%             
%         else
%             
%             peak_angle_diff = NaN;  % if max peak is in corner or center, ignore in analysis
%         end
%         
%         same_or_not_7_plusminus_degree(count2)= peak_angle_diff;
%         
%         count2=count2+1;
%         
%     end
%     
    orientations(h)= new_angle;
    orientations_adjusted(h)= abs(new_angle);
    
end

%figure; hist(same_or_not_7_degree(:));
%figure; hist(same_or_not_7_plusminus_degree(:));


figure; hist(orientations(:));
figure; hist(orientations_adjusted(:));
disp('');