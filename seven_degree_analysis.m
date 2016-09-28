parms.dir_load_data = 'C:\Users\Dori\Desktop\bin size 6 nonsmooth\results updated';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for i=1:length(file_names)
    file_name = file_names{i};
    load(file_name);
        

    %finds the angle between all the points and vertical line from center
    
    six_orientation_pts =S.six_orientation_pts;
    
  m_x= six_orientation_pts(7,1);
    m_y= six_orientation_pts(7,2);
    
    for h= 1:length(six_orientation_pts);
        
        pt_x= six_orientation_pts(h,1);
        pt_y= six_orientation_pts(h,2);
        
        orientation_angles_left(h) = mod((atan2(pt_y-m_y, pt_x-m_x) * 180/pi), 180);
    
        orientation_angles_top(h) = mod(90 - orientation_angles_left(h), 180);
        
        orientation_angles_right(h) = mod(180 - orientation_angles_left(h), 180);
        
        orientation_angles_bottom(h) = mod(180 - orientation_angles_top(h),180);
    end
   
    
    orientation_angles_left(7) = [];
    orientation_angles_right(7) = [];
    orientation_angles_top(7) = [];
    orientation_angles_bottom(7) = [];
  
    
    orientation_left_right=union(orientation_angles_left, orientation_angles_right);
    orientation_top_bottom=union(orientation_angles_top, orientation_angles_bottom);
    
    orientations_all = union(orientation_left_right, orientation_top_bottom);
    
    degree= 7.5;
    tmp= abs(orientation_left_right - degree);
    tmp2=abs(orientation_top_bottom - degree);
    [idx idx] = min(tmp);
    [idex idex]= min(tmp2);
    closest=orientation_left_right(idx);
    closest2=orientation_top_bottom(idex);
    
    temp(1)= abs(closest-degree);
    temp(2)=abs(closest2-degree);
    [idx idx]= min(temp);
    seven_degrees= temp(idx);
    
    all_degrees(i) = seven_degrees;
    
end

figure;
hist(all_degrees(:));
    
figure; hist(orientations_all(:));
    
disp('');





