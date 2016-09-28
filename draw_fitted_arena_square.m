function [draw_fitted_arena] = draw_fitted_arena_square(

parms.dir_load_data1 = 'N:\users\rebekkah\final data smoothed\results with fitted box arena';
%parms.dir_load_data2 = 'N:\users\rebekkah\final data smoothed\major orientation results';

%parms.save_pictures = 'N:\users\rebekkah\final data smoothed\images of fitted box arena';

dir_name= parms.dir_load_data1;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count =1;
count2=1;

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
    
    angle_changes(h)= arena_angle_changes(1);
    
    new_angle= Cell2.smallest_degree + arena_angle_changes(1);
    
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
    
    x_1_r=rotated_sq_x(1);
    x_2_r=rotated_sq_x(2);
    x_3_r=rotated_sq_x(3);
    x_4_r=rotated_sq_x(4);
    
    y_1_r=rotated_sq_y(1);
    y_2_r=rotated_sq_y(2);
    y_3_r=rotated_sq_y(3);
    y_4_r=rotated_sq_y(4);
    
    Cell= Cell2;
    
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    strTitle = sprintf('Cell_i=%d_r%s_d%s_s%s_t%d_c%d',h,Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell);
    fig =figure('name',strTitle);
    % title(strTitle);
    % plot the rat's path (and the spikes on it)
    plot(pos_mean_x,pos_mean_y,'k');hold on;
    %plot(spk_x,spk_y,'.r');
    axis equal;axis off;
    axis ij
    title(file_name);
    hold on;
    
    line([x_1_r, x_3_r],[y_1_r, y_3_r])
    line([x_3_r, x_4_r],[y_3_r, y_4_r])
    line([x_4_r, x_2_r], [y_4_r, y_2_r])
    line([x_2_r, x_1_r], [y_2_r, y_1_r])
end