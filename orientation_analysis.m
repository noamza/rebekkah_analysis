
function [draw_orientation_line, smallest_degree, smallest_degree_wall] = orientation_analysis(rate_mat_after_auto_pearson, rate_mat, six_orientation_pts, draw_orientation_line, max_inds) 
    

    imagesc(rate_mat);axis('xy');axis equal;axis off;hold on;   
    
    axis ij; 

    lw =3;
    
    [size_x, size_y] = size(rate_mat_after_auto_pearson);

    %finds the angle between all the points and vertical line from center
    
    m_x= six_orientation_pts(7,1);
    m_y= six_orientation_pts(7,2);
    
    for h= 1:length(six_orientation_pts);
        
        pt_x= six_orientation_pts(h,1);
        pt_y= six_orientation_pts(h,2);
        
         orientations(h) = atan2(pt_y-m_y, pt_x-m_x) * 180/pi;
        
%         orientation_angles_left(h) = mod((atan2(pt_y-m_y, pt_x-m_x) * 180/pi), 180);
%     
%         orientation_angles_top(h) = mod(90 - orientation_angles_left(h), 180);
%         
%         orientation_angles_right(h) = mod(180 - orientation_angles_left(h), 180);
%         
%         orientation_angles_bottom(h) = mod(180 - orientation_angles_top(h),180);
    end
   
    % lowest angle is 90, removes negative degrees          
   orientations_mod_90= mod(orientations, 90); 
   % lowest angle regardless of location of point
   orientations_abs= 45 - abs(orientations_mod_90-45); 
    %smallest_angle   
    orientations_abs= orientations_abs(1:6);
   smallest_degree= min(orientations_abs);
    
   main_axis_index= find(orientations_abs==smallest_degree);
   main_axis_index= main_axis_index(1); 
   
  original_angle= orientations(main_axis_index);
    
%     orientation_angles_left(7) = [];
%     orientation_angles_right(7) = [];
%     orientation_angles_top(7) = [];
%     orientation_angles_bottom(7) = [];
%     
%     min_top = min(orientation_angles_top);
%     min_bottom = min(orientation_angles_bottom);
%     min_right = min(orientation_angles_right);
%     min_left = min(orientation_angles_left);
%   
    
    [size_x, size_y] = size(rate_mat);
      

    %central points
    x= 1:size_x;
    x2= 1:size_x;
    
    max_inds_len=length(max_inds);
    
    division_arena_x= int64(size_x/3) ;
    division_arena_y= int64(size_y/3) ;
    
    %..........................
    
    % finds index that falls in top-center box when arena divided into 9
    % squares (to use for top line)
    
    inds_y_top= find(max_inds(:,1) > 1 & max_inds(:,1) < division_arena_x);
    inds_y_top=inds_y_top';
    inds_x_top= find(max_inds(:,2) > division_arena_y & max_inds(:,2) < 2*division_arena_y);
    inds_x_top=inds_x_top';
    inds_top= intersect(inds_x_top, inds_y_top);
    
    if ~isempty(inds_top)
        use_pt_top_x= max_inds(inds_top(1),1);
        use_pt_top_y= max_inds(inds_top(1),2);
    elseif ~isempty(inds_y_top)
        use_pt_top_x= max_inds(inds_y_top(1),1);
        use_pt_top_y= max_inds(inds_y_top(1),2); 
    elseif ~isempty(inds_x_top)
        use_pt_top_x= max_inds(inds_x_top(1),1);
        use_pt_top_y= max_inds(inds_x_top(1),2); 
    else
         use_pt_top_x= 0;  % remove error messages
        use_pt_top_y= 0; % remove error messages
    end
    
    
    
   % finds index that falls in bottom-center box when arena divided into 9
    % squares (to use for bottom line)
    
    inds_y_bottom= find(max_inds(:,1) > 2*division_arena_x  & max_inds(:,1) < 3*division_arena_x);
    inds_y_bottom=inds_y_bottom';
    inds_x_bottom= find(max_inds(:,2) > division_arena_y & max_inds(:,2) < 2*division_arena_y);
    inds_x_bottom=inds_x_bottom';
    inds_bottom= intersect(inds_x_bottom, inds_y_bottom);
    
    if ~isempty(inds_bottom)
        use_pt_bottom_x= max_inds(inds_bottom(1),1);
        use_pt_bottom_y= max_inds(inds_bottom(1),2);
    elseif ~isempty(inds_y_bottom)
        use_pt_bottom_x= max_inds(inds_y_bottom(1),1);
        use_pt_bottom_y= max_inds(inds_y_bottom(1),2); 
    elseif ~isempty(inds_x_bottom)
        use_pt_bottom_x= max_inds(inds_x_bottom(1),1);
        use_pt_bottom_y= max_inds(inds_x_bottom(1),2); 
    else 
         use_pt_bottom_x= 0; %remove error messages
        use_pt_bottom_y= 0; %remove error messages
    end
    
    % finds index that falls in center-left box when arena divided into 9
    % squares (to use for left line)
    
    inds_x_left= find(max_inds(:,1) > division_arena_x  & max_inds(:,1) < 2*division_arena_x);
    inds_x_left=inds_x_left';
    inds_y_left= find(max_inds(:,2) > 1 & max_inds(:,2) < division_arena_y);
    inds_y_left=inds_y_left';
    inds_left= intersect(inds_x_left, inds_y_left);
    
    
     if ~isempty(inds_left)
        use_pt_left_x= max_inds(inds_left(1),2);
        use_pt_left_y= max_inds(inds_left(1),1);
    elseif ~isempty(inds_y_left)
        use_pt_left_x= max_inds(inds_y_left(1),2);
        use_pt_left_y= max_inds(inds_y_left(1),1); 
     elseif ~isempty(inds_x_left)
            use_pt_left_x= max_inds(inds_x_left(1),2);
            use_pt_left_y= max_inds(inds_x_left(1),1);
     else
         use_pt_left_x= 0;   % remove error messages
         use_pt_left_y =0;  % to remove error messages
    end
    
    
 
    
    % finds index that falls in center-right box when arena divided into 9
    % squares (to use for right line)
    
    inds_x_right= find(max_inds(:,1) > division_arena_x  & max_inds(:,1) < 2*division_arena_x);
    inds_x_right=inds_x_right';
    inds_y_right= find(max_inds(:,2) > 2*division_arena_y & max_inds(:,2) < 3*division_arena_y);
    inds_y_right=inds_y_right';
    inds_right= intersect(inds_x_right, inds_y_right);
    
    
     if ~isempty(inds_right)
        use_pt_right_x= max_inds(inds_right(1),2);
        use_pt_right_y= max_inds(inds_right(1),1);
    elseif ~isempty(inds_y_right)
        use_pt_right_x= max_inds(inds_y_right(1),2);
        use_pt_right_y= max_inds(inds_y_right(1),1); 
     elseif ~isempty(inds_x_right)
        use_pt_right_x= max_inds(inds_x_right(1),2);
        use_pt_right_y= max_inds(inds_x_right(1),1); 
     else
         use_pt_right_x= 0;
        use_pt_right_y= 0; 
    end
    
    %.........................................
    
     
   if original_angle <= 90 & original_angle >=45 | original_angle >= -135 & original_angle <=-90
        smallest_degree_wall = 4;
       
        %uses center-top index
        y_intercept= -(tand(smallest_degree)*(use_pt_top_y)- use_pt_top_x);
        y= tand(smallest_degree)*x + y_intercept;
        
        %y = tand(min_top)*x ; %top bar y-intercept is 0
        
        plot(x,y, 'w','LineWidth',lw); hold on; % plots top line 
        
        
        %uses center-bottom
        
        y_intercept= -(tand(smallest_degree)*(use_pt_bottom_y)- use_pt_bottom_x);
        y2= tand(smallest_degree)*x2 + y_intercept;
        
         %y_intercept = -(tand(min_top)*(size_x)- size_y); % y-intercept given pt size_x,size_y
         %y2 = tand(min_top)*x2 + y_intercept ;
          plot(x2,y2, 'w','LineWidth',lw); hold on; %plots bottom line
        
          strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n top wall minimum angle', smallest_degree, smallest_degree_wall);
            title(strTitle);
          
          
  elseif original_angle <= 135 & original_angle >=90 | original_angle >= -90 & original_angle <=-45
        smallest_degree_wall = 3;
        
        %y = -tand(min_bottom)*x + size_y;    
        
        y_intercept= (tand(smallest_degree)*(use_pt_top_y)+ use_pt_top_x);
        y= -tand(smallest_degree)*x + y_intercept;
        
        plot(x,y, 'w','LineWidth',lw); hold on;  %plots top line 
        
%          y_intercept = tand(min_bottom)*(size_x);
%          y2 = -tand(min_bottom)*x2 + y_intercept ;

        y_intercept= (tand(smallest_degree)*(use_pt_bottom_y)+ use_pt_bottom_x);
        y2= -tand(smallest_degree)*x2 + y_intercept;
        
          plot(x2,y2, 'w','LineWidth',lw); hold on; %plots bottom line
         
            strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n bottom wall minimum angle', smallest_degree, smallest_degree_wall);
            title(strTitle);
            
            
   elseif original_angle >=135 & original_angle <=180 | original_angle >= -45 & original_angle <=0 | original_angle < -180 
        smallest_degree_wall = 2;
       
         
        % y = -tand(min_right)*x + size_y ;
         
         y_intercept= (tand(smallest_degree)*(use_pt_left_y)+ use_pt_left_x);
        y= -tand(smallest_degree)*x + y_intercept;
        
         plot(y,x, 'w','LineWidth',lw); hold on; %plots left line
         
%           y_intercept = tand(min_right)*(size_x);
%          y2 = -tand(min_right)*x2 + y_intercept ;
        
        y_intercept= (tand(smallest_degree)*(use_pt_right_y)+ use_pt_right_x);
        y2= -tand(smallest_degree)*x2 + y_intercept;

          plot(y2,x2, 'w','LineWidth',lw); hold on;
         title('right wall minimum angle');  
         
           strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n right wall minimum angle', smallest_degree, smallest_degree_wall);
            title(strTitle);
         
    elseif original_angle <= 45 & original_angle >=0 | original_angle <=-135 & original_angle >=-180 | original_angle>= 180
        smallest_degree_wall = 1;      
      % y = tand(min_left)*x ;
      
      y_intercept= -tand(smallest_degree)*(use_pt_left_y)+ use_pt_left_x;
        y= tand(smallest_degree)*x + y_intercept;
        plot(y,x, 'w','LineWidth',lw); hold on;
       
%         y_intercept = -(tand(min_left)*(size_x)- size_y);
%          y2 = tand(min_left)*x2 + y_intercept ;

    y_intercept= -tand(smallest_degree)*(use_pt_right_y)+ use_pt_right_x;
        y2= tand(smallest_degree)*x2 + y_intercept;

        plot(y2,x2, 'w','LineWidth',lw); hold on;
    title('left wall minimum angle');  
        
      strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n left wall minimum angle', smallest_degree, smallest_degree_wall);
            title(strTitle);
    
    else 
        disp('wtf')
    end
    
    
%     plot(x1,y2,'r+');
%     hold on
%     ezplot(eq1);
%     plot(sol.x,sol.y,'ro');
    
     
    %axis ij;
    %strTitle= sprintf('smallest degree = %0.3f \n wall =%d', smallest_degree, smallest_degree_wall)
    %title(strTitle);     
    
    disp('');
    





