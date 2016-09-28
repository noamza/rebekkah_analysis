function [rotated_sq_x, rotated_sq_y] = draw_fitted_arena_box(file_name, dir_load_data,smallest_degree)

parms.dir_load_data1 = 'N:\users\rebekkah\final data smoothed\results with fitted box arena';

count =1;
count2=1;

cd(parms.dir_load_data1);
dat1= load(file_name);
cd(dir_load_data);
Cell1=dat1.S;

% for h =1:length(file_names)
%     cd(parms.dir_load_data1);
%     file_name = file_names{h};
%     dat1 = load(file_name);
%     Cell1= dat1.S;
%     
%     cd(parms.dir_load_data2);
%     dat2 = load(file_name);
%     Cell2= dat2.S;
   

    max_pts_inside= find(Cell1.phi_of_arena.percent_inside_final(:,1)==max(Cell1.phi_of_arena.percent_inside_final(:,1)));
    
    max_pts_inside= max_pts_inside(1);
    
    arena_angle_changes= Cell1.phi_of_arena.percent_inside_final(max_pts_inside,5);
        
    %new_angle= smallest_degree + arena_angle_changes(1);
    
    a= Cell1.phi_of_arena.percent_inside_final(max_pts_inside, 4);
    a= sqrt(a);
    x_0= Cell1.phi_of_arena.percent_inside_final(max_pts_inside, 2)/100;
    y_0= Cell1.phi_of_arena.percent_inside_final(max_pts_inside, 3)/100;
    phi= Cell1.phi_of_arena.percent_inside_final(max_pts_inside, 5);
    
    x_min= x_0-a/2;
    x_max= x_0+a/2;
    y_min= y_0-a/2;
    y_max=y_0+a/2;
    
    %%rotating square
    
    square(1:2,1)= [x_min y_min];
    square(1:2,2)= [x_min y_max];
    square(1:2,3)= [x_max y_min];
    square(1:2,4)= [x_max y_max];
    
    square_x_pts= [x_min x_min x_max x_max];
    square_y_pts= [y_min y_max y_min y_max];
    
    rotation_matrix= [cosd(phi) -sind(phi); sind(phi) cosd(phi)];
    
    for j=1:4;
        rotated_square(1:2,j)=rotation_matrix*square(1:2,j);
        rotated_sq_x(j)=rotated_square(1,j);
        rotated_sq_y(j)=rotated_square(2,j);
    end
    
   
            
end