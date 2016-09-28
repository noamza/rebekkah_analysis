
function [min_angle, main_axis, axis_and_orient] = adjusted_orientation_analysis(rate_mat_after_auto_pearson, rate_mat, six_orientation_pts, phi) 
    

  %  imagesc(rate_mat);axis('xy');axis equal;axis off;hold on;   
    
 %   axis ij; 

%    lw =3;
    
    [size_x, size_y] = size(rate_mat_after_auto_pearson);

%     auto_distances = [];
%     new_distances = [];
% 
%     %finds distances from center of auto_corr_map to all max peaks
%     
%     cen = 1:length(auto_max_inds);
%     auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
% 
%     new_distances = sort(auto_distances (:));
%     
%     % finds the 6 closest pts to the center
%     
%     auto_dist_inds1 = [];
%     auto_dist_inds2 = [];
%     auto_dist_inds3 = [];
%     auto_dist_inds4 = [];
%     auto_dist_inds5 = [];
%     auto_dist_inds6 = [];
%    
%     auto_dist_inds1 = find(auto_distances == new_distances(2));
%     auto_dist_inds2 = find(auto_distances == new_distances(3));
%     auto_dist_inds3 = find(auto_distances == new_distances(4));
%     auto_dist_inds4 = find(auto_distances == new_distances(5));
%     auto_dist_inds5 = find(auto_distances == new_distances(6));
%     auto_dist_inds6 = find(auto_distances == new_distances(7));
%     
%     union1 = union(auto_dist_inds1, auto_dist_inds2);
%     union2 = union(auto_dist_inds3, auto_dist_inds4);
%     union3 = union(auto_dist_inds5, auto_dist_inds6);
%     union1 = union(union1, union2);
%     
%     auto_dist_inds = union(union1, union3);
 
    
%     for h= 1:length(auto_dist_inds);
%         six_orientation_pts (h,1) = auto_max_inds(auto_dist_inds(h),1);
%         six_orientation_pts (h,2) = auto_max_inds(auto_dist_inds(h),2);
%     end
%     
    %finds the angle between all the points and vertical line from center
    
    m_x= six_orientation_pts(7,1); % center point x
    m_y= six_orientation_pts(7,2); % center point y
    
    %finds orientations -180 to 180
    
    for h= 1:length(six_orientation_pts); 
        
        pt_x= six_orientation_pts(h,1);
        pt_y= six_orientation_pts(h,2);
        
        orientations(h) = atan2(pt_y-m_y, pt_x-m_x) * 180/pi;
    end

    % accounts arena angle difference 
       orientations_arena=orientations+phi;
    % removes unnessary last point (center point)
       orientations_arena=orientations_arena(1:6);
    % lowest angle is 90, removes negative degrees          
   orientations_mod_90= mod(orientations_arena, 90); 
   % lowest angle regardless of location of point
   orientations_abs= 45 - abs(orientations_mod_90-45); 
    %smallest_angle   
   min_angle= min(orientations_abs);
   
   main_axis_index= find(orientations_abs==min_angle);
   main_axis_index= main_axis_index(1); 
   
  original_angle= orientations_arena(main_axis_index);
  
  %finds axis of smallest angle (north-south or east-west) and also whether
  %negative or positive slope of smallest angle axis
  
  if original_angle <= 45 & original_angle >=0 | original_angle <=-135 & original_angle >-180 | original_angle> 180
      main_axis= 1; %south-north axis
      axis_and_orient= 1;%negative slope to south-north axis
  elseif original_angle <= 90 & original_angle >=45 | original_angle >= -135 & original_angle <=-90
   main_axis= 0; %east-west axis
   axis_and_orient=4; %negative east-west
  elseif original_angle <= 135 & original_angle >=90 | original_angle >= -90 & original_angle <=-45
     main_axis= 0; %east-west axis
   axis_and_orient=3; %positive east-west
  elseif original_angle >=135 & original_angle <=180 | original_angle >= -45 & original_angle <=0 | original_angle <= -180 
     main_axis= 1; %north-south axis
   axis_and_orient=2; %positive south-north
  else
      disp('WTF. Something is wrong.');
  end
   
        disp('')
        
        
%         orientation_angles_left(h) = mod((atan2(pt_y-m_y, pt_x-m_x) * 180/pi), 360);
%     
%         orientation_angles_top(h) = mod(90 - orientation_angles_left(h), 360);
%         
%         orientation_angles_right(h) = mod(180 - orientation_angles_left(h), 360);
%         
%         orientation_angles_bottom(h) = mod(180 - orientation_angles_top(h),360);
    
   
    % calculate points to the right and left, and top and bottom to
    % determine which is cardinal axis and which is positive and negative
    % axis
      
    %find angles to the right of the center
    
%     pts_to_the_bottom_of_center= find(six_orientation_pts(:,1)> m_x);
%     pts_to_the_top_of_center= find(six_orientation_pts(:,1)< m_x);
%     pts_to_the_right_of_center= find(six_orientation_pts(:,2)> m_y);
%     pts_to_the_left_of_center= find(six_orientation_pts(:,2)< m_y);
%     
%     %bottom points correct. orientations not.
%     
%     pt_y_bottom= six_orientation_pts(pts_to_the_bottom_of_center,2); 
%     pt_x_bottom= six_orientation_pts(pts_to_the_bottom_of_center,1); 
%     
%     bottom_orientations= mod((atan2(pt_y_bottom -m_y, pt_x_bottom -m_x) * 180/pi), 360);
%     bottom_orientations= bottom_orientations + phi;
%     
%     % left works accurately 
%     
%     pt_y_left= six_orientation_pts(pts_to_the_left_of_center,2); 
%     pt_x_left= six_orientation_pts(pts_to_the_left_of_center,1); 
%     
%     left_orientations= mod((atan2(pt_y_left -m_y, pt_x_left -m_x) * 180/pi), 180);
%     left_orientations= left_orientations + phi;
%     
%     pt_y_top= six_orientation_pts(pts_to_the_top_of_center,2); 
%     pt_x_top= six_orientation_pts(pts_to_the_top_of_center,1); 
%     
%     pt_y_right= six_orientation_pts(pts_to_the_right_of_center,2); 
%     pt_x_right= six_orientation_pts(pts_to_the_right_of_center,1); 
%     
%    right_orientations= mod((90- atan2(pt_y_right -m_y, pt_x_right -m_x) * 180/pi), 180);
%    right_orientations= right_orientations + phi;
%     
%     %% see if difference between 180 and 360
%     
%     for h= 1:length(six_orientation_pts);
%         
%         pt_x= six_orientation_pts(h,1);
%         pt_y= six_orientation_pts(h,2);
%         
%         orient_180_angles_left(h) = mod((atan2(pt_y-m_y, pt_x-m_x) * 180/pi), 180);
%     
%         orient_180_angles_top(h) = mod(90 - orient_180_angles_left(h), 180);
%         
%         orient_180_angles_right(h) = mod(180 - orient_180_angles_left(h), 180);
%         
%         orient_180_angles_bottom(h) = mod(180 - orient_180_angles_top(h),180);
%     end
%     
%     
%     
%     orientation_angles_left_adj= orientation_angles_left + phi;
%     orientation_angles_left_adj= mod(orientation_angles_left_adj, 360)
%     orientation_angles_top_adj= orientation_angles_top - phi;
%     orientation_angles_top_adj= mod(orientation_angles_top_adj, 360)
%     orientation_angles_right_adj= orientation_angles_right - phi;
%     orientation_angles_right_adj= mod(orientation_angles_right_adj, 360)
%     orientation_angles_bottom_adj= orientation_angles_bottom+ phi;
%     orientation_angles_bottom_adj= mod(orientation_angles_bottom_adj, 360)
%     
%     orient_180_angles_left= orient_180_angles_left + phi;
%     orient_180_angles_left= mod(orient_180_angles_left, 180)
%     orient_180_angles_top= orient_180_angles_top - phi;
%     orient_180_angles_top= mod(orient_180_angles_top, 180)
%     orient_180_angles_right= orient_180_angles_right - phi;
%     orient_180_angles_right= mod(orient_180_angles_right, 180)
%     orient_180_angles_bottom= orient_180_angles_bottom+ phi;
%     orient_180_angles_bottom= mod(orient_180_angles_bottom, 180)
%     
%     
%     orientation_angles_left(7) = [];
%     orientation_angles_right(7) = [];
%     orientation_angles_top(7) = [];
%     orientation_angles_bottom(7) = [];
%     
%     orientation_angles_left_adj(7) = [];
%     orientation_angles_right_adj(7) = [];
%     orientation_angles_top_adj(7) = [];
%     orientation_angles_bottom_adj(7) = [];
%     
%     orient_180_angles_left(7)= [];
%     orient_180_angles_right(7)= []
%     orient_180_angles_top(7)=[];
%     orient_180_angles_bottom(7)=[];
%     
%     min_top = min(orientation_angles_top);
%     min_bottom = min(orientation_angles_bottom);
%     min_right = min(orientation_angles_right);
%     min_left = min(orientation_angles_left);
%   
%     min_top_adj = min(orientation_angles_top_adj);
%     min_bottom_adj = min(orientation_angles_bottom_adj);
%     min_right_adj = min(orientation_angles_right_adj);
%     min_left_adj = min(orientation_angles_left_adj);
%   
%     min_top_180 = min(orient_180_angles_top);
%     min_bottom_180 = min(orient_180_angles_bottom);
%     min_right_180 = min(orient_180_angles_right);
%     min_left_180 = min(orient_180_angles_left);
%     
%     [size_x, size_y] = size(rate_mat);
%       
% 
%     %central points
%     x= 1:size_x;
%     x2= 1:size_x;
%     
%     max_inds_len=length(max_inds);
%     
%     division_arena_x= int64(size_x/3) ;
%     division_arena_y= int64(size_y/3) ;
%     
%     %..........................
%     
%     % finds index that falls in top-center box when arena divided into 9
%     % squares (to use for top line)
%     
%     inds_y_top= find(max_inds(:,1) > 1 & max_inds(:,1) < division_arena_x);
%     inds_y_top=inds_y_top';
%     inds_x_top= find(max_inds(:,2) > division_arena_y & max_inds(:,2) < 2*division_arena_y);
%     inds_x_top=inds_x_top';
%     inds_top= intersect(inds_x_top, inds_y_top);
%     
%     if ~isempty(inds_top)
%         use_pt_top_x= max_inds(inds_top(1),1);
%         use_pt_top_y= max_inds(inds_top(1),2);
%     elseif ~isempty(inds_y_top)
%         use_pt_top_x= max_inds(inds_y_top(1),1);
%         use_pt_top_y= max_inds(inds_y_top(1),2); 
%     else
%         use_pt_top_x= max_inds(inds_x_top(1),1);
%         use_pt_top_y= max_inds(inds_x_top(1),2); 
%     end
%     
%     
%     
%    % finds index that falls in bottom-center box when arena divided into 9
%     % squares (to use for bottom line)
%     
%     inds_y_bottom= find(max_inds(:,1) > 2*division_arena_x  & max_inds(:,1) < 3*division_arena_x);
%     inds_y_bottom=inds_y_bottom';
%     inds_x_bottom= find(max_inds(:,2) > division_arena_y & max_inds(:,2) < 2*division_arena_y);
%     inds_x_bottom=inds_x_bottom';
%     inds_bottom= intersect(inds_x_bottom, inds_y_bottom);
%     
%     if ~isempty(inds_bottom)
%         use_pt_bottom_x= max_inds(inds_bottom(1),1);
%         use_pt_bottom_y= max_inds(inds_bottom(1),2);
%     elseif ~isempty(inds_y_bottom)
%         use_pt_bottom_x= max_inds(inds_y_bottom(1),1);
%         use_pt_bottom_y= max_inds(inds_y_bottom(1),2); 
%     else
%         use_pt_bottom_x= max_inds(inds_x_bottom(1),1);
%         use_pt_bottom_y= max_inds(inds_x_bottom(1),2); 
%     end
%     
%     % finds index that falls in center-left box when arena divided into 9
%     % squares (to use for left line)
%     
%     inds_x_left= find(max_inds(:,1) > division_arena_x  & max_inds(:,1) < 2*division_arena_x);
%     inds_x_left=inds_x_left';
%     inds_y_left= find(max_inds(:,2) > 1 & max_inds(:,2) < division_arena_y);
%     inds_y_left=inds_y_left';
%     inds_left= intersect(inds_x_left, inds_y_left);
%     
%     
%      if ~isempty(inds_left)
%         use_pt_left_x= max_inds(inds_left(1),2);
%         use_pt_left_y= max_inds(inds_left(1),1);
%     elseif ~isempty(inds_y_left)
%         use_pt_left_x= max_inds(inds_y_left(1),2);
%         use_pt_left_y= max_inds(inds_y_left(1),1); 
%     else
%             use_pt_left_x= max_inds(inds_x_left(1),2);
%             use_pt_left_y= max_inds(inds_x_left(1),1); 
%     end
%     
%     
%  
%     
%     % finds index that falls in center-right box when arena divided into 9
%     % squares (to use for right line)
%     
%     inds_x_right= find(max_inds(:,1) > division_arena_x  & max_inds(:,1) < 2*division_arena_x);
%     inds_x_right=inds_x_right';
%     inds_y_right= find(max_inds(:,2) > 2*division_arena_y & max_inds(:,2) < 3*division_arena_y);
%     inds_y_right=inds_y_right';
%     inds_right= intersect(inds_x_right, inds_y_right);
%     
%     
%      if ~isempty(inds_right)
%         use_pt_right_x= max_inds(inds_right(1),2);
%         use_pt_right_y= max_inds(inds_right(1),1);
%     elseif ~isempty(inds_y_right)
%         use_pt_right_x= max_inds(inds_y_right(1),2);
%         use_pt_right_y= max_inds(inds_y_right(1),1); 
%      else
%         use_pt_right_x= max_inds(inds_x_right(1),2);
%         use_pt_right_y= max_inds(inds_x_right(1),1); 
%     end
%     
%     %.........................................
%     
%      
%     if min_top <= min_right & min_top <= min_left & min_top <= min_bottom
%         smallest_degree_wall = 1;
%         smallest_angle = min_top;
%         
%         %uses center-top index
%         y_intercept= -(tand(min_top)*(use_pt_top_y)- use_pt_top_x);
%         y= tand(min_top)*x + y_intercept;
%         
%         %y = tand(min_top)*x ; %top bar y-intercept is 0
%         
%         plot(x,y, 'w','LineWidth',lw); hold on; % plots top line 
%         
%         
%         %uses center-bottom
%         
%         y_intercept= -(tand(min_top)*(use_pt_bottom_y)- use_pt_bottom_x);
%         y2= tand(min_top)*x2 + y_intercept;
%         
%          %y_intercept = -(tand(min_top)*(size_x)- size_y); % y-intercept given pt size_x,size_y
%          %y2 = tand(min_top)*x2 + y_intercept ;
%           plot(x2,y2, 'w','LineWidth',lw); hold on; %plots bottom line
%         
%           strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n top wall minimum angle', smallest_angle, smallest_degree_wall);
%             title(strTitle);
%           
%           
%     elseif min_bottom <= min_right & min_bottom <= min_left & min_bottom <= min_bottom
%         smallest_degree_wall = 3;
%         smallest_angle = min_bottom;
%         
%         %y = -tand(min_bottom)*x + size_y;    
%         
%         y_intercept= (tand(min_bottom)*(use_pt_top_y)+ use_pt_top_x);
%         y= -tand(min_bottom)*x + y_intercept;
%         
%         plot(x,y, 'w','LineWidth',lw); hold on;  %plots top line 
%         
% %          y_intercept = tand(min_bottom)*(size_x);
% %          y2 = -tand(min_bottom)*x2 + y_intercept ;
% 
%         y_intercept= (tand(min_bottom)*(use_pt_bottom_y)+ use_pt_bottom_x);
%         y2= -tand(min_bottom)*x2 + y_intercept;
%         
%           plot(x2,y2, 'w','LineWidth',lw); hold on; %plots bottom line
%          
%             strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n bottom wall minimum angle', smallest_angle, smallest_degree_wall);
%             title(strTitle);
%             
%             
%     elseif min_right <= min_top & min_right <= min_left & min_right <= min_bottom 
%         smallest_degree_wall = 4;
%          smallest_angle = min_right;
%          
%         % y = -tand(min_right)*x + size_y ;
%          
%          y_intercept= (tand(min_right)*(use_pt_left_y)+ use_pt_left_x);
%         y= -tand(min_right)*x + y_intercept;
%         
%          plot(y,x, 'w','LineWidth',lw); hold on; %plots left line
%          
% %           y_intercept = tand(min_right)*(size_x);
% %          y2 = -tand(min_right)*x2 + y_intercept ;
%         
%         y_intercept= (tand(min_right)*(use_pt_right_y)+ use_pt_right_x);
%         y2= -tand(min_right)*x2 + y_intercept;
% 
%           plot(y2,x2, 'w','LineWidth',lw); hold on;
%          title('right wall minimum angle');  
%          
%            strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n right wall minimum angle', smallest_angle, smallest_degree_wall);
%             title(strTitle);
%          
%     elseif min_left <= min_top & min_left <= min_right & min_left <= min_bottom 
%         smallest_degree_wall = 2;
%         smallest_angle = min_left;
%         
%       % y = tand(min_left)*x ;
%       
%       y_intercept= -tand(min_left)*(use_pt_left_y)+ use_pt_left_x;
%         y= tand(min_left)*x + y_intercept;
%         plot(y,x, 'w','LineWidth',lw); hold on;
%        
% %         y_intercept = -(tand(min_left)*(size_x)- size_y);
% %          y2 = tand(min_left)*x2 + y_intercept ;
% 
%     y_intercept= -tand(min_left)*(use_pt_right_y)+ use_pt_right_x;
%         y2= tand(min_left)*x2 + y_intercept;
% 
%         plot(y2,x2, 'w','LineWidth',lw); hold on;
%     title('left wall minimum angle');  
%         
%       strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n left wall minimum angle', smallest_angle, smallest_degree_wall);
%             title(strTitle);
%     
%     else 
%         disp('wtf')
%     end
%     
    
%     plot(x1,y2,'r+');
%     hold on
%     ezplot(eq1);
%     plot(sol.x,sol.y,'ro');
    
     
    %axis ij;
    %strTitle= sprintf('smallest degree = %0.3f \n wall =%d', smallest_degree, smallest_degree_wall)
    %title(strTitle);
    
   %...............find smallest adjusted for arena angle difference grid orientation angle....
   
%     if min_top_adj <= min_right_adj & min_top_adj <= min_left_adj & min_top_adj <= min_bottom_adj
%         smallest_degree_adj = min_top_adj;
%         smallest_degree_wall_adj = 1;
%     elseif min_bottom_adj <= min_right_adj & min_bottom_adj <= min_left_adj & min_bottom_adj <= min_top_adj
%         smallest_degree_adj = min_bottom_adj;
%         smallest_degree_wall_adj = 3;
%     elseif min_left_adj <= min_right_adj & min_left_adj <= min_top_adj & min_left_adj <= min_bottom_adj
%         smallest_degree_adj = min_left_adj;
%         smallest_degree_wall_adj= 2;
%     elseif min_right_adj <= min_top_adj & min_right_adj <= min_left_adj & min_right_adj <= min_bottom_adj
%         smallest_degree_adj= min_right_adj;
%         smallest_degree_wall_adj= 4;
%     end
%  %..........checks to see if taking 180 or 360 cutoffs makes a differnce- [it doesnt, as it shouldnt]... 
%     
%        if min_top_180 <= min_right_180 & min_top_180 <= min_left_180 & min_top_180 <= min_bottom_180
%         smallest_degree_180 = min_top_180;
%         smallest_degree_wall_180 = 1;
%     elseif min_bottom_180 <= min_right_180 & min_bottom_180 <= min_left_180 & min_bottom_180 <= min_top_180
%         smallest_degree_180 = min_bottom_180;
%         smallest_degree_wall_180 = 3;
%     elseif min_left_180 <= min_right_180 & min_left_180 <= min_top_180 & min_left_180 <= min_bottom_180
%         smallest_degree_180 = min_left_180;
%         smallest_degree_wall_180= 2;
%     elseif min_right_180 <= min_top_180 & min_right_180 <= min_left_180 & min_right_180 <= min_bottom_180
%         smallest_degree_180= min_right_180;
%         smallest_degree_wall_180= 4;
%     end
    
    disp('');
    





