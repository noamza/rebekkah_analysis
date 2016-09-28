function FindMinStretchShiftRescaling 

parms.dir_load_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\max';
parms.dir_save_data = 'C:\Users\Dori\Desktop\rescalign arena results for shuffling';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

 rates_first= cell(1,47);
 rates_new= cell(1,47);

% enumerate on cells
for j =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{j};
    load(file_name);
    
    %   load('214_06-03-31.mat_cell1.mat')
    
    max_inds1_orig= FindMaxIndsRateMap(rate_mats.arena1);
    autocorr1= Cross_Correlation(rate_mats.arena1,rate_mats.arena1);
    autoMaxInds1= FindAutoMaxInds(autocorr1);
    PF_radius1= findPlaceFieldRadius(autocorr1, autoMaxInds1);
    max_inds1_orig = RemoveTooCloseMaxIndsRescaling(max_inds1_orig, PF_radius1, rate_mats.arena1);
    zone_mat_orig= CreateZoneMat(rate_mats.arena1, rate_mats.arena1, max_inds1_orig);
    
    false_field_thresh= nanmean2(rate_mats.arena1);
    zone_mat_orig(zone_mat_orig<= false_field_thresh)= 0;
    
    peak_rates = findPeakRates(max_inds1_orig, zone_mat_orig);
    false_ind= find(peak_rates<= false_field_thresh);
    max_inds1_orig(false_ind, :)=[];
    
    max_inds1_orig_orig= max_inds1_orig;
    
    for arena_count= 2:4
        
        max_inds1_orig= max_inds1_orig_orig;
        
        
        name= sprintf('arena%d', arena_count);
        
        max_inds2_orig= FindMaxIndsRateMap(rate_mats.(name));
        autocorr2= Cross_Correlation(rate_mats.(name),rate_mats.(name));
        autoMaxInds2= FindAutoMaxInds(autocorr2);
        PF_radius2= findPlaceFieldRadius(autocorr2, autoMaxInds2);
        max_inds2_orig = RemoveTooCloseMaxIndsRescaling(max_inds2_orig, PF_radius2, rate_mats.(name));
        zone_mat_change= CreateZoneMat(rate_mats.(name), rate_mats.(name), max_inds2_orig);
        
        false_field_thresh= nanmean2(rate_mats.(name));
        zone_mat_change(zone_mat_change<= false_field_thresh)= 0;
        
        peak_rates = findPeakRates(max_inds2_orig, zone_mat_change);
        false_ind= find(peak_rates<= false_field_thresh);
        max_inds2_orig(false_ind, :)=[];
        
        
        %testing below%
        %    max_inds1_orig= [1 1; 30 30; 15 15;];
        %   max_inds2_orig= [1 0.5; 15 7.5; 30 15];
        
        %testing above%
        
        max_inds1= max_inds1_orig;
        max_inds2=max_inds2_orig;
        
        size_1= length(rate_mats.arena1);
        size_2= length(rate_mats.arena1);
        
        size_max= max(size_1,size_2);
        %try all As from 1 to 2
        %try all Bs from 0 to size of arena
        
        %size_max=50;
        
        [len1 ~]= size(max_inds1);
        [len2 ~]= size(max_inds2);
        
        min_len= min(len1, len2);
        
        
        len_diff= len1/len2;
        
        
        % if abs(1- (1/len_diff)) < 0.3  & len2~= 1
        
        if len1 > 2 && len2 > 2
            
            
            %for pairs= 1:min_len
            
            % change to smaller --> larger arena always
            % make matrix of cells vs arena config by stretch (color) -->
            % normalize each cell by mean stretch
            
            % M: ls. vr. hr. ss.
            
            min_dist=inf;
            %  rates(h)= new_max_inds2;
            % zone_vals(h) = zone_val;
            % count=count+1;
            
            A_stretches= -1.2;
           
            
            for A=2.^(-1:0.2:1) 
                
                A_stretches= A_stretches+0.2;
                
                for B=0:10:size_max  
                    
                     C_stretches= -1.2;
                    
                    for C=2.^(-1:0.2:1)  
                        
                         C_stretches= C_stretches+0.2;
                        
                        for D=0:10:size_max
                            dist=0;
                            
                            x1= max_inds2_orig(:,1);
                            y1= max_inds2_orig(:,2);
                            
                            for k=1:length(x1)
                                max_inds2(k,1)= A*(x1(k))-B;
                                max_inds2(k,2)= C*(y1(k))-D;
                            end
                            %% finds closest pairs
                            
                            max_inds1= max_inds1_orig;
                            
                            
                            zone_val= nan(1,len1);
                            
                            
                            
                            for h=1:min_len
                                
                                [min_i, min_j] = FindTransfPtToPt2(max_inds1, max_inds2);
                                
                                new_max_inds2(min_j,1)= max_inds1(min_i,1);
                                new_max_inds2(min_j,2)= max_inds1(min_i,2);
                                
                                max_inds1(min_i,:)= nan;
                                max_inds2(min_j,:)= nan;
                                
                                zone_val(min_i)= min_j;
                            end
                            
                            nan_inds= find(isnan(zone_val));
                            
                            try
                                max_inds1_orig(nan_inds, :)= [];
                            end
                            
                            zone_val(isnan(zone_val))= [];
                            
                            x= max_inds1_orig(:,1);
                            y=max_inds1_orig(:,2);
                            
                            for i=1:length(x)
                                
                                max_inds2= max_inds2_orig(zone_val, :);
                                
                                x1= max_inds2(:,1);
                                y1= max_inds2(:,2);
                                
                                %%how to decide minimum distance
                                
                                dist= dist + (x(i) - A*x1(i) -B)^2 + (y(i) - C*y1(i) -D)^2;
                            end
                            
                            dist= dist/length(x);
                            
                            if dist< min_dist;
                                min_dist= dist;
                                A_min=A;
                                B_min=B;
                                C_min=C;
                                D_min=D;
                                
                                A_s_m= A_stretches;
                                C_s_m= C_stretches;
                                
                                zone_val_min=zone_val;
                                
                                %  n_min=n;
                                %  m_min=m;
                                
                            end
                            
                        end
                    end
                end
            end
            
            
            %         end
            %end
            
            
            %     max_inds1= max_inds1_orig(perms_n(n_min,:),:);
            %    max_inds2= max_inds2_orig(perms_m(m_min,:),:);
            
            if min_dist <= 70
                
                x= max_inds1_orig(:,1);
                y= max_inds1_orig(:,2);
                
                x1= max_inds2_orig(:,1);
                y1= max_inds2_orig(:,2);
                
                final_x1= nan(1,min_len);
                final_y1= nan(1,min_len);
                for i=1:min_len
                    final_x1(i)= A_min*(x1(i))+B_min;
                    final_y1(i)= C_min*(y1(i))+D_min;
                end
                
                max_inds_final= nan(min_len, 2);
                max_inds_final(:,1)= final_x1;
                max_inds_final(:,2)= final_y1;
                
                max_inds2_final= nan(min_len,2);
                
                %    figure;
                %    plot(max_inds1_orig(:,1), max_inds1_orig(:,2), 'o', 'MarkerFaceColor', 'k'); hold on;
                %    plot(max_inds_final(:,1), max_inds_final(:,2), 'o', 'MarkerFaceColor', 'r');
                
                [sizex sizey] = size(rate_mats.arena1);
                [size_x size_y]= size(rate_mats.(name));
                
                x_diff= sizex/size_x;
                y_diff= sizey/size_y;
                
                stretch_x(count)= A_min/x_diff;
                shift_x(count)=B_min;
                stretch_y(count)= C_min/y_diff;
                shift_y(count)=D_min;
                
                A_all(count)= A_s_m;
                B_all(count)= B_min;
                C_all(count)= C_s_m;
                D_all(count)= D_min;
                
                
                
                [orig_x orig_y]= size(zone_mat_orig);
                [new_x new_y]= size(zone_mat_change);
                
                factorx= A_min;
                factory= C_min;
                
                tform = maketform('affine',[factory 0 0; 0 factorx 0; 0 0 1]);
                stretch_arena = imtransform(zone_mat_change,tform);
                
                
                max_inds_final= round(max_inds_final);
                
                % [sorted_final zone_val_min]= sort(zone_val_min);
                
                % max_inds_final(1:min_len,[1 2])= max_inds_final(zone_val_min,[1 2]);
                
                max_inds2_final(1:min_len,[1 2])= max_inds2_orig(zone_val_min,[1 2]);
                
                %plot images
                
                
                
                
                % to plot figures
%                 
%                                 figure;
%                                 subplot(1,2,1);
%                                 imagesc(zone_mat_orig); hold on;
%                 
%                                 color_val= [{'bx'}, {'kx'}, {'rx'}, {'gx'}, {'mx'}, {'wx'}, {'yx'}, {'bo'}, {'ko'}, {'ro'}, {'go'}, {'mo'} , {'wo'}] ;
%                 
%                                 for cen=1:min_len
%                                     plot(max_inds1_orig(cen,2),max_inds1_orig(cen,1), color_val{cen}, 'MarkerSize', 10,'LineWidth', 3);
%                                 end
%                 
%                                 title(sprintf('%0.2f', min_dist));
%                 
%                                 subplot(1,2,2);
%                                 imagesc(zone_mat_change); hold on;
%                                 for cen=1:min_len
%                                     plot(max_inds2_final(cen,2),max_inds2_final(cen,1), color_val{cen}, 'MarkerSize', 10,'LineWidth', 3);
%                                 end
                
                [max_inds_len ~]= size(max_inds_final);
                
                rates_orig= nan(1,max_inds_len);
                rates_change= nan(1,max_inds_len);
                
                for h= 1:max_inds_len;
                    rates_orig(h)= zone_mat_orig(max_inds1_orig(h,1),max_inds1_orig(h,2));
                    rates_change(h)= zone_mat_change(max_inds2_final(h,1),max_inds2_final(h,2));
                end
                
                corrcoef_rate= corrcoef(rates_orig, rates_change);
                
              
                
                %                 cd(parms.dir_save_data);
                %                 file= sprintf('%s_arena%d.mat', file_name, arena_count);
                %                 save(file, 'rates_orig', 'rates_change');
                %                 cd(parms.dir_load_data);
                
%                 figure;
%                 n=1;
%                 m=3;
%                 
%                 subplot(n,m,1)
%                 imagesc(rate_mats.arena1);
%                 title(sprintf('%d %d', sizex, sizey));
%                 subplot(n,m,2)
%                 imagesc(rate_mats.(name)); 
%                 title(sprintf('%d %d', size_x, size_y));
%                 
%                 subplot(n,m,3);
%  
%                 plot(max_inds1_orig_orig(:,2),max_inds1_orig_orig(:,1), 'x', 'MarkerSize', 10,'LineWidth', 3);
%                 hold on;
%                 
%                 title(sprintf('%0.2f %0.2f %0.2f %0.2f', ...
%                     stretch_x(count), stretch_y(count), A_all(count), C_all(count)));
%       
%              
%                     plot(max_inds2_final(:,2),max_inds2_final(:,1), 'rx', 'MarkerSize', 10,'LineWidth', 3);
%             
%                 axis ij
                
               
                
                 if length(corrcoef_rate) > 1
                    
                    corr_rates(count)= corrcoef_rate(2);
                    
                    rates_first{count}= rates_orig;
                    rates_new{count}= rates_change;
                    
                    count=count+1;
             
                end
                                
                
            end
            
        end
        
    end
    
    
    
    
    disp('');
    
    clearvars -except file_names j parms ...
        stretch_x shift_x stretch_y shift_y count corr_rates A_all C_all rates_first rates_new
    
end

save('corr rates of fitted arenas', 'corr_rates', 'rates_first', 'rates_new')


figure; hist(stretch_x);
figure; hist(stretch_y);

figure; hist(A_all);
figure; hist(C_all);

figure; hist(corr_rates);

mean(corr_rates)

disp('');
