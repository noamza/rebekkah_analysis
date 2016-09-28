
cd('N:\users\rebekkah\results and info of analysis')
load('rescaling arenas info')

count=1;

for i= 1:length(all_rate_mats)  % open info for every cell
    
    %original arena info
    max_indices= all_max_inds{i};
    max_inds_orig= max_indices{1}; %max inds of first arena
    peak_rates_all_arena= all_peak_rates{i};
    peak_rates_orig= peak_rates_all_arena{1};
    rate_mats_all_arena= all_rate_mats{i};
    rate_mat_orig= rate_mats_all_arena{1};
    zone_mats_all_arena= all_zone_mats{i};
    
    num_of_arenas= length(peak_rates_all_arena);
    
    [len1 ~]= size(max_inds_orig); %num of PF in original arena
    
    max_inds_orig_orig= max_inds_orig;
        
    for arena_count= 2:num_of_arenas-1 %for all except last arena
        
        [len2 ~]= size(max_indices{arena_count}); %num of PF in next arena
        len_diff= len1/len2;    %ratio diff of PF nums
        
        
        
        %if diff isnt too large and if more than one PF in next arena
        if abs(1- (1/len_diff)) < 0.3  & len2~= 1
            
            max_inds_orig= max_inds_orig_orig; %reset to remove nans from minimum dist analysis 
       
            rate_mat_next= rate_mats_all_arena{arena_count};
            
            size_1= length(rate_mat_orig);
            size_2= length(rate_mat_next);
            size_max= max(size_1,size_2); % finds larger arena size for shift
                   
            max_inds_next= max_indices{arena_count};
            
            %testing BELOW
%             max_inds_orig=[0 0; 20 20; 10 10; 26 26];
%             max_inds_next=[52 52; 20 20; 40 40; 0 0; 9 10];
%             
%             size_max=50;
%             
%             len1= 4; 
%             len2= 5;
            %testing ABOVE
            
            min_dist=inf;
            
            min_len= min(len1, len2);
            
            x_norm_stretch= -1.2; %x stretch
            for A=2.^(-1:0.2:1)
                x_norm_stretch= x_norm_stretch+0.2;
                
                for B=0:10:size_max  %x shift
                    
                    y_norm_stretch= -1.2;
                    for C=2.^(-1:0.2:1)  %y stretch
                        y_norm_stretch= y_norm_stretch+0.2;

                        for D=0:10:size_max % y shift
                            
                            x1= max_inds_next(:,1);    %x of next arena
                            y1= max_inds_next(:,2);    %y of next arena
                             
                            %applies stretch and shift to next arena coords
                            max_inds2= nan(size(max_inds_next));
                            for k=1:length(x1)
                                max_inds2(k,1)= A*(x1(k))-B;    %add stretch and scale
                                max_inds2(k,2)= C*(y1(k))-D;
                            end
                            
                            % finds closest pair of points
                            max_inds1= max_inds_orig; 
                            
                            zone_val= nan(1,len2); % the number zonemat in order to orig arena pt pairs
                            
                            new_max_inds2=nan(min_len,2); 
                            for h=1:min_len
                                
                                %two pt pairs with most min distance between
                                [min_i, min_j] = FindTransfPtToPt2(max_inds1, max_inds2);
                                
                                new_max_inds2(min_i,:)= max_inds2(min_j,:); 
                                
                                max_inds1(min_i,:)= nan;
                                max_inds2(min_j,:)= nan;
                                
                                zone_val(min_j)= min_i; %the order of 2nd inds so it'll pair with 1st inds
                            end
                            
                            zone_missing= setdiff(1:len1, zone_val);
                            
                            try
                                max_inds_orig(zone_missing, :)= []; %removes unnecessary index if there is
                            end
                            
                            % zone_val(isnan(zone_val))= [];
                            
                            x= max_inds_orig(:,1);
                            y= max_inds_orig(:,2);
                            
                             x1= new_max_inds2(:,1);
                             y1= new_max_inds2(:,2);
                                
                             dist=0;
                            for h=1:length(x)
                                %finds distance and then compares to min
                                %distance
                                dist= dist+ Distance(x(h),y(h),x1(h),y1(h));
                            end
                            
                            dist= dist/length(x);
                            
                            if dist< min_dist;
                                min_dist= dist;
                                A_min=A;
                                B_min=B;
                                C_min=C;
                                D_min=D;
                                
                                min_x_st_norm= x_norm_stretch;
                                min_y_st_norm= y_norm_stretch;
                                
                                zone_val_min=zone_val;

                                %check if this is same as final inds
                                min_max_inds2= new_max_inds2;
                                
                                % zone_val= nan(1,len1);
                            end
                        end
                    end
                end    
            end
            
            
              %testing to make sure same as min_max_inds2:correct 
%                 x= max_inds_orig(:,1);
%                 y= max_inds_orig(:,2);
%                 
%                 x1= max_inds2_orig(:,1);
%                 y1= max_inds2_orig(:,2);
%                             
%                 final_x1= nan(1,min_len);
%                 final_y1= nan(1,min_len);
%                for h=1:min_len
%                     final_x1(h)= A_min*(x1(h))-B_min;
%                     final_y1(h)= C_min*(y1(h))-D_min;
%                 end 
%                   
%                 max_inds_final= nan(min_len, 2);
%                 max_inds_final(:,1)= final_x1;
%                 max_inds_final(:,2)= final_y1;
%                 max_inds_final= max_inds_final(zone_val_min,:); %put in correct order 
                
                % figure
                % plot rate maps, zone maps, and indices
                fig=figure;
                n=2;
                m=4;
                subplot(n,m,1)
                imagesc(rate_mat_orig);
                subplot(n,m,5)
                imagesc(rate_mat_next);
                subplot(n,m,2)
                imagesc(zone_mats_all_arena{1});
                subplot(n,m,6)
                imagesc(zone_mats_all_arena{arena_count});
            
                subplot(n,m,[3 4 7 8]);
                plot(max_inds_orig(:,1), max_inds_orig(:,2), 'ko', 'MarkerFaceColor', 'k'); hold on;
                plot(min_max_inds2(:,1), min_max_inds2(:,2), 'ro', 'MarkerFaceColor', 'r');
                for h=1:length(max_inds_orig)
                plot([max_inds_orig(h,1);min_max_inds2(h,1)], ...
                    [max_inds_orig(h,2);min_max_inds2(h,2)], '-k');
                end
                axis ij;
                
                
                cd('N:\users\rebekkah\results and info of analysis\rescaling best fit images')
                saveas(fig,sprintf('recaling num %d.jpg', count));
                
                all_x_stretches(count)= min_x_st_norm;  
                all_y_stretches(count)= min_y_st_norm;
                
                count=count+1;
                
                close all
                
                disp('');
        end
        
    end

    disp('')
    
end

disp('');