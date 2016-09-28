
dbstop if error

cd('N:\users\rebekkah\results and info of analysis')
load('rescaling arenas info.mat')

load('file names rescaling data.mat')

count=1;

for i=1:length(all_rate_mats)
    
    max_indices= all_max_inds{i};
    rate_mats= all_rate_mats{i};
    peak_rates_all_arenas= all_peak_rates{i};
    zone_mats=all_zone_mats{i};
    
      max_inds_orig= max_indices{1}; %
    [orig_x orig_y]= size(rate_mats{1}); %first arena dimensions
    
    [len1,~]= size(max_inds_orig);
    
    for arena_count= 2:length(max_indices)-1 %opens all but last arena
        
        max_inds_orig= max_indices{1}; % max inds of first arena
        rate_mat_orig= rate_mats{1};
        
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
%             tform = maketform('affine',[factory 0 0; 0 factorx 0; 0 0 1]);
%             stretch_arena = imtransform(zone_mats{arena_count},tform);
            
            max_inds_stretch= max_inds_next;
            max_inds_stretch(:,1)= max_inds_next(:,1) * factorx;
            max_inds_stretch(:,2)= max_inds_next(:,2) * factory;
            
            min_len= min(len1,len2); % number of least fields
            
            zone_val= nan([1 len1]);
            
            new_max_inds_next=nan([len2 2]);
            
            %find pairs of field centers closest together
            for h=1:min_len
                
                [min_i, min_j] = FindTransfPtToPt2(max_inds_orig, max_inds_stretch);

                max_inds_orig(min_i,:)= nan;
                max_inds_stretch(min_j,:)= nan;
                
                zone_val(min_j)= min_i;
            end
            
            peak_rates_next=nan(1,len1); 
            rates=peak_rates_all_arenas{arena_count};
            for h=1:length(zone_val);
                if ~isnan(zone_val(h))
            peak_rates_next(zone_val(h))=rates(h);
                end
            end
            
            zone_missing= setdiff(1:len1, zone_val);
         
            peak_rates_orig= peak_rates_all_arenas{1};
            
            try 
            peak_rates_orig(zone_missing)=nan;
            end
            
         rates_orig{count}= peak_rates_orig;
         rates_next{count}=peak_rates_next;
            
            % TO REMOVE MAX FIELD uncomment below:
%                          max_ind= find(peak_rates_orig==max(peak_rates_orig));
%                          peak_rates_orig(max_ind)=[];
%                          peak_rates_next(max_ind)=[];
%             
%                           if length(peak_rates_orig) > 1
%                          rates_corr= corrcoef(peak_rates_orig, peak_rates_next);
%                          rates_corr=rates_corr(2);
%                          else
%                              rates_corr=nan;
%                          end
            % TO REMOVE MAX FIELD uncomment above:
            
            nan_ind= find(isnan(peak_rates_orig));
            peak_rates_orig(nan_ind)=[];
            peak_rates_next(nan_ind)=[];
            
             rates_corr= corrcoef(peak_rates_orig, peak_rates_next);
             rates_corr=rates_corr(2);
            
            title= sprintf('%s arena %d', filenames{i}, arena_count);
            filename_rescaling{count}= title;

            rates_corrcoef(count)= rates_corr;
                       
            count=count+1;
            
        end
        
    end
end

save('filenames of rescaling firing comparison', 'filename_rescaling')

 % save('firing stability rescaled arenas wo hyperfield', 'rates_orig', 'rates_next', 'rates_corrcoef')

%     length(rates_corr_all)
%
%     rates_corr_mean(k)= nanmean(rates_corr_all);
%     rates_corr_sum(k)= sum(rates_corr_all >= 0.5);
%     rates_corr_sum2(k)= sum(rates_corr_all >= 0.7);
%
%   %  clear rates_corr_all;
% end
%
% save('shuffled_remapping_rate_stability first and last arena wo max2', 'rates_corr_mean', 'rates_corr_sum', 'rates_corr_sum2')
%
% %real results
% % mean =0.4609
% % sum > 0.5 = 61
% % sum > 0.6 = 58

