%%% look at trajectory to see how to fit square (for walls)

parms.dir_load_data = 'C:\Users\Dori\Desktop\testing top lines';
parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\data sets\fitted arena with HD without circular rooms_ 202 cells';


dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    dat= load(file_name);
    Cell=dat.S;
    
    % calculate the the rat's head direction (using average of both leds)
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
%             strTitle = sprintf('Cell_i=%d_r%s_d%s_s%s_t%d_c%d',h,Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell);
%                 fig =figure('name',strTitle);
%                % title(strTitle);
%                 % plot the rat's path (and the spikes on it)
%                 plot(pos_mean_x,pos_mean_y,'k');hold on;
%                 %plot(spk_x,spk_y,'.r');
%                 axis equal;axis off;
%                 axis ij
%                 title(file_name);
%                 hold on;
    
    x_0_est_median= median(pos_mean_x);
    y_0_est_median= median(pos_mean_y);
    
    sorted_x= sort(pos_mean_x);
    x_len= length(pos_mean_x);
    max_x_percentile= int64(0.999*x_len);
    min_x_percentile= int64(0.001*x_len);
    
    x_max_est= sorted_x(max_x_percentile);
    x_min_est= sorted_x(min_x_percentile);
    
    x_line = x_max_est-x_min_est;
    x_0_est= x_min_est + x_line/2;
    
    sorted_y= sort(pos_mean_y);
    y_len= length(pos_mean_y);
    max_y_percentile= int64(0.999*y_len);
    min_y_percentile= int64(0.001*y_len);
    
    y_max_est= sorted_y(max_y_percentile) ;
    y_min_est= sorted_y(min_y_percentile);
    
    y_line = y_max_est-y_min_est;
    y_0_est= y_min_est+  y_line/2;
    
    x_y_a_est(1)= x_max_est-x_min_est;
    x_y_a_est(2)= y_max_est-y_min_est;
    
    a_est=mean(x_y_a_est);
    
    x= x_min_est:x_max_est;
    y= y_min_est:y_max_est;
    
    line([x_min_est, x_min_est],[y_min_est, y_max_est])
    line([x_max_est, x_max_est],[y_min_est, y_max_est])
    line([x_min_est, x_max_est], [y_min_est, y_min_est])
    line([x_min_est, x_max_est], [y_max_est, y_max_est])
    
    plot(x_0_est, y_0_est, 'x', 'markersize', 15, 'color', 'r');
    
    % x_0_est, y_0_est
    % a_est
    % phi
    
    
    %%
    
    xv = [x_min_est x_min_est x_max_est x_max_est x_min_est];
    yv = [y_max_est y_min_est y_min_est y_max_est y_max_est];
    
    xq= pos_mean_x;
    yq= pos_mean_y;
    
    inside= inpolygon(xq,yq,xv,yv);
    sum_inside= sum(inside);
    pos_len= length(pos_mean_x);
    
    percent_inside_est = sum_inside/pos_len;
    
    count=1;
    
    disp('');
    
    for x_0= x_0_est-5: x_0_est+5;
        for y_0= y_0_est-5: y_0_est+5;
            for a= 0.85*a_est : 1.15*a_est;
                for phi= -5:5
                    
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
                    
                    
                    %side= x_max_est-x_min_est;
                    
                    %  y_diff = side*tand(abs(phi));
                    
                    
                    xv = [x_1_r x_3_r x_4_r x_2_r x_1_r];
                    yv = [y_1_r y_3_r y_4_r y_2_r y_1_r];
                    
                    xq= pos_mean_x;
                    yq= pos_mean_y;
                    
                    inside= inpolygon(xq,yq,xv,yv);
                    sum_inside= sum(inside);
                    pos_len= length(pos_mean_x);
                    
                    percent_inside = sum_inside/pos_len;
                    
                    area= a^2;
                    
                    percent_inside_x_y_a_phi(count,1)= percent_inside*100;
                    percent_inside_x_y_a_phi(count,2)= x_0*100;
                    percent_inside_x_y_a_phi(count,3)= y_0*100;
                    percent_inside_x_y_a_phi(count,4)= area/100;
                    percent_inside_x_y_a_phi(count,5)= phi;
                    
                    count=count+1;
                    
                    disp('');
                    
                end
                
                disp('');
            end
            disp('');
            
        end
        
    end
    
    disp('')
    
    less_than_99= find(percent_inside_x_y_a_phi(:,1) < 99);
    
    percent_inside_x_y_a_phi(less_than_99,:) = [];
    
    min_area = find(percent_inside_x_y_a_phi(:,4) == min(percent_inside_x_y_a_phi(:,4)))
    min_area=min_area'
    
    percent_inside_final= percent_inside_x_y_a_phi(min_area,:);
    
    
%       strTitle = sprintf('Cell_i=%d_r%s_d%s_s%s_t%d_c%d',h,Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell);
%                 fig =figure('name',strTitle);
%                % title(strTitle);
%                 % plot the rat's path (and the spikes on it)
%                 plot(pos_mean_x,pos_mean_y,'k');hold on;
%                 %plot(spk_x,spk_y,'.r');
%                 axis equal;axis off;
%                 axis ij
%                 title(file_name);
%                 hold on;
%     
%         line([x_1_r, x_3_r],[y_1_r, y_3_r])
%         line([x_3_r, x_4_r],[y_3_r, y_4_r])
%         line([x_4_r, x_2_r], [y_4_r, y_2_r])
%         line([x_2_r, x_1_r], [y_2_r, y_1_r])
                
                
    disp('');
    
    S.phi_of_arena.percent_inside_final= percent_inside_final;
    Cell.phi_of_arena.percent_inside_final=S.phi_of_arena.percent_inside_final;
    
    %saves results
    cd(parms.dir_save_data);
    title1 =strcat (file_name);
    save (title1,'S');
    cd(parms.dir_load_data);
    
end

