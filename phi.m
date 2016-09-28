function phi

parms.dir_load_data = 'C:\Users\Dori\Desktop\final data smoothed\results';
parms.dir_save_pictures= 'C:\Users\Dori\Desktop\final data smoothed\major orientation';
parms.dir_save_data = 'C:\Users\Dori\Desktop\final data smoothed\major orientation results'; 


dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

lw =3;

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
    six_orientation_pts=S.six_orientation_pts;
    
    six_orientation_pts(7,:)=[];
    
    count=1;
    
   for points= 1:length(six_orientation_pts);
       for points2= points+1:length(six_orientation_pts)-1
        distance(count,1)= Distance(six_orientation_pts(points,1), six_orientation_pts(points,2),...
       six_orientation_pts(points2,1), six_orientation_pts(points2,2));
        distance(count,2)= points;
        distance(count,3)= points2;
        count=count+1;
       end
   end
    
   [max_value, index] = max(distance(:,1));
   pt1= distance(index,2); 
   pt2= distance(index,3);
      
   m_x= six_orientation_pts(pt1,1);
   m_y= six_orientation_pts(pt1,2);
            
        pt_x= six_orientation_pts(pt2,1);
        pt_y= six_orientation_pts(pt2,2);
        
        orientation_angle_left = mod((atan2(pt_y-m_y, pt_x-m_x) * 180/pi), 180);
        orientation_angle_top= mod(90 - orientation_angle_left, 180); 
        orientation_angle_right = mod(180 - orientation_angle_left, 180);
        orientation_angle_bottom = mod(180 - orientation_angle_top,180);
        
        fig=figure;
        subplot(1,3,2)
        imagesc(S.rate_mat);axis equal;axis off;hold on;
        
        %central points
        [size_x, size_y] =size(S.rate_mat);
        
        x= 1:size_x;
        x2= 1:size_x;
        
        if orientation_angle_left <= orientation_angle_right & orientation_angle_left<=orientation_angle_top & orientation_angle_left<= orientation_angle_bottom
            angle=orientation_angle_left;
            wall = 2;
            
           y = tand(angle)*x ;
            plot(y,x, 'w','LineWidth',lw); hold on;
       
            y_intercept = -(tand(angle)*(size_x)- size_y);
            y2 = tand(angle)*x2 + y_intercept ;
            plot(y2,x2, 'w','LineWidth',lw); hold on; 
            
            
        elseif orientation_angle_right <= orientation_angle_left & orientation_angle_right<=orientation_angle_top & orientation_angle_right<=orientation_angle_bottom
            angle=orientation_angle_right;
            wall=4;
            
            y = -tand(angle)*x + size_y ;
         
            plot(y,x, 'w','LineWidth',lw); hold on;
         
            y_intercept = tand(angle)*(size_x);
             y2 = -tand(angle)*x2 + y_intercept ;
            plot(y2,x2, 'w','LineWidth',lw); hold on;  
            
            
        elseif orientation_angle_bottom<=orientation_angle_top 
            angle= orientation_angle_bottom;
            wall= 3;
            
              y = -tand(angle)*x + size_y;
        
        plot(x,y, 'w','LineWidth',lw); hold on;
        
         y_intercept = tand(angle)*(size_x);
         y2 = -tand(angle)*x2 + y_intercept ;
          plot(x2,y2, 'w','LineWidth',lw); hold on; 
          
        elseif orientation_angle_top<=orientation_angle_bottom
            angle= orientation_angle_top;
            wall=1;
            
           y = tand(angle)*x ;
        
        plot(x,y, 'w','LineWidth',lw); hold on;
        
         y_intercept = -(tand(angle)*(size_x)- size_y);
         y2 = tand(angle)*x2 + y_intercept ;
          plot(x2,y2, 'w','LineWidth',lw); hold on;
            
                   
           
        else
            disp('huh')
        end
        
         
       
        fig_pair= subplot(1,3,1);
            fig_pair= plot_ellipse(S.autocorr, S.module.major, S.module.minor ,S.module.phi, S.six_orientation_pts, fig_pair);

     
            subplot(1,3,3);
            imagesc(S.zone_mat); axis equal;
        
        angles(h)= angle;
        angles_left(h) =orientation_angle_left;
        angles_top(h)=orientation_angle_top;
        
     
    cd(parms.dir_save_pictures);
       saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',S.rat,S.date,S.session,S.tetrode,S.cell)); %         % debugger - return
        saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',S.rat,S.date,S.session,S.tetrode,S.cell)); % 
cd(parms.dir_load_data);
        
close all        

S.major_angle=angle;
S.wall=wall; 

cd(parms.dir_save_data);  
title1 =strcat (file_name);
save (title1,'S');
cd(parms.dir_load_data);
        
    disp('');
        
end

figure; hist(angles);


