

parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';

% parms.dir_load_data2= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\cluster matrixs above 0.3';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count_all=0;


 count1=0;
 count2=0;
 
% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    A=load(file_name);
    %    cd(parms.dir_load_data2);
    %     file_name2= sprintf('%s.mat', file_name);
    %     B= load(file_name2);
    
    
    
    rate_mat = A.rate_mats.arena1;
    
    autocorr= Cross_Correlation(rate_mat, rate_mat); %creates autocorrelation of rate mat
    
    PF_radius= findPlaceFieldRadius(autocorr);   %finds the radius of the field
    
    six_orientation_pts = find_six_points(autocorr, PF_radius); %finds the six points of the hexagon
    
    size_rate_mat= size(rate_mat);
    
    x_arena_len= size_rate_mat(2);
    y_arena_len= size_rate_mat(1);
    horizon_len= max(six_orientation_pts(:,2)) - min(six_orientation_pts(:,2));
    vertical_len= max(six_orientation_pts(:,1)) - min(six_orientation_pts(:,1));
    
%     x_arena_len_diff=nan(1,4);
%     y_arena_len_diff=nan(1,4);
%     horizon_len_diff=nan(1,4);
%     vertical_len_diff=nan(1,4);
    
    arena1_size=size(rate_mat);
    
    count=1;
    
    six_orientation_pts=[];
    
    if isfield(A.rate_mats,'arena5')
        
        for arena_count=2:5;
            
            fldnm = sprintf('arena%d', arena_count);
            rate_mat = A.rate_mats.(fldnm);
            
            autocorr= Cross_Correlation(rate_mat, rate_mat); %creates autocorrelation of rate mat
            
            PF_radius= findPlaceFieldRadius(autocorr);   %finds the radius of the field
            
            six_orientation_pts = find_six_points(autocorr, PF_radius); %finds the six points of the hexagon
            
            size_rate_mat= size(rate_mat);
            
            x_arena_len_diff(1)= size_rate_mat(2)/arena1_size(2); %percentage of xaxis change 
            y_arena_len_diff(1)= size_rate_mat(1)/arena1_size(1); %percentage of yaxis change
            horizon_len_diff(1)=  (max(six_orientation_pts(:,2)) - min(six_orientation_pts(:,2)))/horizon_len;
            vertical_len_diff(1)= (max(six_orientation_pts(:,1)) - min(six_orientation_pts(:,1)))/vertical_len; 
            
            count=count+1;
            
            count_all=count_all+1;
           
            
             if x_arena_len_diff> 0.9 & x_arena_len_diff < 1.1       %if x-axis is the same
                count1=count1+1;
                unchanged_axis(count1)= horizon_len_diff/x_arena_len_diff;   
                unchanged_other_axis(count1)= vertical_len_diff/x_arena_len_diff;
          elseif x_arena_len_diff< 0.9 | x_arena_len_diff > 1.1         %if x-axis is different
              count2=count2+1;
                manipulated_axis(count2)= horizon_len_diff/x_arena_len_diff;
                manipulated_other_axis(count2)= vertical_len_diff/x_arena_len_diff;
            end
          
        if y_arena_len_diff> 0.9 & y_arena_len_diff < 1.1
                count1=count1+1;
                unchanged_axis(count1)= vertical_len_diff/y_arena_len_diff;   
                unchanged_other_axis(count1)= horizon_len_diff/y_arena_len_diff;
          elseif y_arena_len_diff< 0.9 | y_arena_len_diff > 1.1
              count2=count2+1;
                 manipulated_axis(count2)= vertical_len_diff/y_arena_len_diff; 
                 manipulated_other_axis(count2)= horizon_len_diff/y_arena_len_diff;
        end
              
                
            crl_x_horiz(count_all)= horizon_len_diff/x_arena_len_diff;
            crl_y_horiz(count_all)= horizon_len_diff/y_arena_len_diff;
            crl_x_vert(count_all)= vertical_len_diff/x_arena_len_diff;
            crl_y_vert(count_all)= vertical_len_diff/y_arena_len_diff;
            
            
        end
        
    else
        
        for arena_count=2:4;
            
            fldnm = sprintf('arena%d', arena_count);
            rate_mat = A.rate_mats.(fldnm);
            
            autocorr= Cross_Correlation(rate_mat, rate_mat); %creates autocorrelation of rate mat
            
            PF_radius= findPlaceFieldRadius(autocorr);   %finds the radius of the field
            
            six_orientation_pts = find_six_points(autocorr, PF_radius); %finds the six points of the hexagon
            
            size_rate_mat= size(rate_mat);
                        
            x_arena_len_diff(1)= size_rate_mat(2)/arena1_size(2); %percentage of xaxis change 
            y_arena_len_diff(1)= size_rate_mat(1)/arena1_size(1); %percentage of yaxis change
            horizon_len_diff(1)=  (max(six_orientation_pts(:,2)) - min(six_orientation_pts(:,2)))/horizon_len;
            vertical_len_diff(1)= (max(six_orientation_pts(:,1)) - min(six_orientation_pts(:,1)))/vertical_len; 
            
            count=count+1;
            
            count_all=count_all+1;
           
            
            if x_arena_len_diff> 0.9 & x_arena_len_diff < 1.1       %if x-axis is the same
                count1=count1+1;
                unchanged_axis(count1)= horizon_len_diff/x_arena_len_diff;   
                unchanged_other_axis(count1)= vertical_len_diff/x_arena_len_diff;
          elseif x_arena_len_diff< 0.9 | x_arena_len_diff > 1.1         %if x-axis is different
              count2=count2+1;
                manipulated_axis(count2)= horizon_len_diff/x_arena_len_diff;
                manipulated_other_axis(count2)= vertical_len_diff/x_arena_len_diff;
            end
          
        if y_arena_len_diff> 0.9 & y_arena_len_diff < 1.1
                count1=count1+1;
                unchanged_axis(count1)= vertical_len_diff/y_arena_len_diff;   
                unchanged_other_axis(count1)= horizon_len_diff/y_arena_len_diff;
          elseif y_arena_len_diff< 0.9 | y_arena_len_diff > 1.1
              count2=count2+1;
                 manipulated_axis(count2)= vertical_len_diff/y_arena_len_diff; 
                 manipulated_other_axis(count2)= horizon_len_diff/y_arena_len_diff;
        end
              
                
            crl_x_horiz(count_all)= horizon_len_diff/x_arena_len_diff;
            crl_y_horiz(count_all)= horizon_len_diff/y_arena_len_diff;
            crl_x_vert(count_all)= vertical_len_diff/x_arena_len_diff;
            crl_y_vert(count_all)= vertical_len_diff/y_arena_len_diff;
            
            
        end
        
        
        
    end
    
    %   figure; plot(x_arena_len_diff, horizon_len_diff, 'x', 'MarkerSize', 15);
    %   title('length of x arena to horizontal spacing');
    %   lsline;
    %  figure; plot(y_arena_len_diff, horizon_len_diff, 'x','MarkerSize', 15);
    %  title('length of y arena to horizontal spacing')
    %  lsline;
    %  figure; plot(x_arena_len_diff, vertical_len_diff, 'x', 'MarkerSize', 15);
    %  lsline;
    %  title('lenth of x arena to vertical spacing')
    %   figure; plot(y_arena_len_diff, vertical_len_diff, 'x', 'MarkerSize', 15);
    %   lsline;
    %   title('length of y arena to vertical spacing')
    
    %     count_all=count_all+1;
    %
    %     %correl_len = corrcoef(x_arena_len_diff, horizon_len_diff);
    %     crl_x_horiz(count_all)= horizon_len_diff/x_arena_len_diff;
    %     %correl_len = corrcoef(y_arena_len_diff, horizon_len_diff);
    %     crl_y_horiz(count_all)= horizon_len_diff/y_arena_len_diff;
    %     %correl_len = corrcoef(x_arena_len_diff, vertical_len_diff);
    %     crl_x_vert(count_all)= vertical_len_diff/x_arena_len_diff;
    %     %correl_len = corrcoef(y_arena_len_diff, vertical_len_diff);
    %     crl_y_vert(count_all)= vertical_len_diff/y_arena_len_diff;
    
    
end

save('grid spacing vs arena deformation above 1', 'crl_x_horiz', 'crl_y_horiz', 'crl_x_vert', 'crl_y_vert', 'manipulated_axis', 'unchanged_axis', ...
            'manipulated_other_axis', 'unchanged_other_axis');

 plot(1,mean(manipulated_axis),'*m'); hold on;
plot(3,mean(unchanged_axis),'*g') ;
plot(2,mean(manipulated_other_axis), '*r')
plot(4,mean(unchanged_other_axis), '*k')
errorbar(1,mean(manipulated_axis),std(manipulated_axis),'.m','linewidth',1);
errorbar(3,mean(unchanged_axis),std(unchanged_axis), '.g','linewidth',1);
errorbar(2,mean(manipulated_other_axis),std(manipulated_other_axis),'.m','linewidth',1);
errorbar(4,mean(unchanged_other_axis),std(unchanged_other_axis), '.g','linewidth',1);
