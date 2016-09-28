
dbstop if error

cd('N:\users\rebekkah\results and info of analysis')
load('rescaling arenas info.mat')

count=1;

for i=1:length(all_rate_mats)
    
    %properties of the cell
    max_indices= all_max_inds{i};
    peak_rates_all_arenas= all_peak_rates{i};
    zone_mats=all_zone_mats{i};
    rate_mats=all_rate_mats{i};
    
    %properities of the original arena
    [orig_x orig_y]= size(rate_mats{1}); %first arena dimensions
    max_inds_orig= max_indices{1}; % max inds of first arena
    rate_mat_orig= rate_mats{1};
   
    [len1,~]= size(max_inds_orig);
    
   % for arena_count= 2:length(max_indices)-1 %opens all but last arena
   for arena_count= length(max_indices) %opens last arena 
   
         peak_rates_orig=peak_rates_all_arenas{1}; %should be reset after each arena_count loop
        max_inds_next= max_indices{arena_count};
        
        [len2,~]= size(max_inds_next);
        len_diff= len1/len2;    %ratio diff of PF nums
        
        %if diff isnt too large and if more than one PF in next arena
        % len2~=1 for regular analysis and >2 for hyperfield excluded
        % analysis
        if abs(1- (1/len_diff)) < 0.3  && len2~= 1
            
            [new_x new_y]= size(zone_mats{arena_count});
            
            factorx= orig_x/new_x;
            factory= orig_y/new_y;
            
            %rescales zone mat to same size as original arena
             tform = maketform('affine',[factory 0 0; 0 factorx 0; 0 0 1]);
             stretch_arena = imtransform(zone_mats{arena_count},tform);
            
             %to image stretch_arena
%              figure;
%              subplot(1,3,1)
%              imagesc(zone_mats{1});
%             % title('size= %d %d', orig_x, orig_y);
%              subplot(1,3,2)
%              imagesc(zone_mats{arena_count});
%             % title('size= %d %d', size(zone_mats{arena_count}));
%              subplot(1,3,3)
%              imagesc(stretch_arena);
%             % title('size= %d %d', size(stretch_arena));
             
            peak_rates_next=nan(1,length(max_inds_orig));
             for h=1:length(max_inds_orig);
                 peak_rates_next(h)=stretch_arena(max_inds_orig(h,1),max_inds_orig(h,2));
             end
             
             zero_inds= find((peak_rates_next==0));
             peak_rates_next(zero_inds)=[];
             peak_rates_orig(zero_inds)=[];
                   
         rates_orig{count}= peak_rates_orig;
         rates_next{count}=peak_rates_next;
            
            % TO REMOVE MAX FIELD uncomment below:
                         max_ind= find(peak_rates_orig==max(peak_rates_orig));
                         peak_rates_orig(max_ind)=[];
                         peak_rates_next(max_ind)=[];
            
                          if length(peak_rates_orig) > 1
                         rates_corr= corrcoef(peak_rates_orig, peak_rates_next);
                         rates_corr=rates_corr(2);
                         else
                             rates_corr=nan;
                         end
            % TO REMOVE MAX FIELD uncomment above:
            
          
            
          %   rates_corr= corrcoef(peak_rates_orig, peak_rates_next);
           %  rates_corr=rates_corr(2);
            
            rates_corrcoef(count)= rates_corr;
            count=count+1;
            
            
            
        end
        
    end
end

save('firing stability same arenas wo hyperfield', 'rates_orig', 'rates_next', 'rates_corrcoef')


