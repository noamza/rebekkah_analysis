function [med_fired, mean_fired, med_nonfired, mean_nonfired] = FindClosestPtCOM(start_pass, end_pass, place_field_t, pos_mean_x, pos_mean_y, field_zero_inds, zero_inds, non_zero_inds, max_inds, PF_radius)
    
%PF_radius_conv= 2.7* PF_radius;

size= max(pos_mean_x) - min(pos_mean_x);
bin_size=3;

conv_max_inds= max_inds*bin_size - (size/2); %convert from rate_mat coords to plot coords

max_inds= conv_max_inds;

for h= 1:length(start_pass)
    
    
    com_y= max_inds(place_field_t(start_pass(h)),1);
    com_x= max_inds(place_field_t(start_pass(h)),2);

    
    count=1;
    distances=[];
    
    for j= start_pass(h): end_pass(h)
    
        
        if ~any(field_zero_inds==j)
        distances(count)= Distance(pos_mean_x(j), pos_mean_y(j), com_x, com_y); 
        
        count=count+1;
        
        elseif any(field_zero_inds==j)
            distances(count) = nan;
            
            count=count+1;
            
        end
    end
    
    min_ind= find(distances== min(distances));
    min_ind=min_ind-1;
    min_ind=min_ind(1); 
    
    pt_x(h)= pos_mean_x(start_pass(h)+min_ind);
    pt_y(h)= pos_mean_y(start_pass(h)+min_ind);
    
    min_distances(h)= min(distances);
    
end

closest_pt_fired_x= pt_x;
closest_pt_fired_y= pt_y;
closest_pt_fired_x(zero_inds)= [];
closest_pt_fired_y(zero_inds)= [];
distances_fired= min_distances;
distances_fired(zero_inds)=[];

closest_pt_non_fired_x= pt_x;
closest_pt_non_fired_y= pt_y;
closest_pt_non_fired_x(non_zero_inds)= [];
closest_pt_non_fired_y(non_zero_inds)= [];
distances_nonfired= min_distances;
distances_nonfired(non_zero_inds)=[];


mean_fired= mean(distances_fired);
med_fired= median(distances_fired);

mean_nonfired= mean(distances_nonfired);
med_nonfired= median(distances_nonfired);

% n=2;
% m=4;
% 
% subplot(n,m,5)
%  plot(closest_pt_non_fired_x,closest_pt_non_fired_y, '.');hold on;
%  
% % plot circles
%  
%   for h=1:length(max_inds)
%         circle(max_inds(h,2),max_inds(h,1), PF_radius_conv)
%   end
%  
%     axis ij; axis equal; axis off;
%     title('closest pass point of no firing')
%     
% 
%     
% 
% subplot(n,m,6)
%  plot(closest_pt_fired_x,closest_pt_fired_y, '.');hold on; 
%  
%  for h=1:length(max_inds)
%         circle(max_inds(h,2),max_inds(h,1), PF_radius_conv)
%   end
%  
%     axis ij; axis equal; axis off;
%     title('closest pass point of firing')