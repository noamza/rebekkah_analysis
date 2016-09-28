
%%open files 

parms.dir_load_data2 = 'N:\users\rebekkah\final data smoothed\results with fitted box arena';
parms.dir_load_data1 = 'N:\users\rebekkah\final data smoothed\FINAL ADJUSTED ANGLES';

dir_name= parms.dir_load_data1;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};


for i =1:length(file_names)
    cd(parms.dir_load_data1);
    file_name = file_names{i};
    dat1 = load(file_name);
    Cell1= dat1.S;
    
    cd(parms.dir_load_data2);
    dat2 = load(file_name);
    Cell2= dat2.S;
    
   
% find smallest degree among 4 boundaries 

    six_orientation_pts = find_six_points(Cell1.autocorr, Cell1.PF_radius)

    six_orientation_pts
    
    m_x= six_orientation_pts(7,1);
    m_y= six_orientation_pts(7,2);
    
    for h= 1:length(six_orientation_pts);
        
        pt_x= six_orientation_pts(h,1);
        pt_y= six_orientation_pts(h,2);
        
        orientation_angles_left(h) = mod((atan2(pt_y-m_y, pt_x-m_x) *180/pi), 180);
    
        orientation_angles_top(h) = mod(90 - orientation_angles_left(h), 180);
        
        orientation_angles_right(h) = mod(180 - orientation_angles_left(h), 180);
        
        orientation_angles_bottom(h) = mod(180 - orientation_angles_top(h),180);
    end
   
    % takes into account arena disorientation 
    
   max_pts_inside= find(Cell2.phi_of_arena.percent_inside_final(:,1)==max(Cell2.phi_of_arena.percent_inside_final(:,1)));
 
    
   phi = Cell2.phi_of_arena.percent_inside_final(max_pts_inside,5);
   phi= phi(1);
    
    orientation_angles_left_adj= orientation_angles_left - phi;
    orientation_angles_left_adj= mod(orientation_angles_left_adj, 180)
    orientation_angles_top_adj= orientation_angles_top + phi;
    orientation_angles_top_adj= mod(orientation_angles_top_adj, 180)
    orientation_angles_right_adj= orientation_angles_right + phi;
    orientation_angles_right_adj= mod(orientation_angles_right_adj, 180)
    orientation_angles_bottom_adj= orientation_angles_bottom- phi;
    orientation_angles_bottom_adj= mod(orientation_angles_bottom_adj, 180) 

    all_orientations(1,1:6)= orientation_angles_left_adj(1:6);
    all_orientations(2,1:6)= orientation_angles_top_adj(1:6);
    all_orientations(3,1:6)= orientation_angles_right_adj(1:6);
    all_orientations(4,1:6)= orientation_angles_bottom_adj(1:6);
    
    [orientation_of_6_pts,indices] = min(all_orientations);
    
% see if there is correlation between max peak location and offset of
% points 
%(0 and 30 should have no max peak-- no sheering)

check_sheering(i, 1:6) = orientation_of_6_pts;


end

save('check sheering 180', 'check_sheering'); 


