load('norm ranks rescaling.mat')
load('rescaling arenas info.mat')

% for k=1:1000
% 
 count=1;
for i=1:length(norm_ranks)

    norm_ranks_arena= norm_ranks{i};  %opens the norm ranks for this cell
    peak_rates_arena=all_peak_rates{i};
    norm_inds_arena= all_norm_inds{i};
       
   % norm_ranks_arena(end)= []; %remove last arena info
    sum_less= sum(norm_ranks_arena <= 0.6); % counts number of arena where hyperfield is not more than half of the most spike count 
    
    if sum_less >= 2 % if hyperfield not firing more than half in at least 2 of the arenas
        
        for arena_count=1:length(peak_rates_arena)
            norm_inds=norm_inds_arena{arena_count};
            max_ind= find(max(peak_rates_arena{arena_count}));
            hyperfield_inds(arena_count,:)= norm_inds(max_ind,:);
            
            i
        end 
            
        inds= find(norm_ranks_arena<=0.6);
        
        hyper_pt_distances= nan(length(inds));
        
       field_max_inds= hyperfield_inds(inds,:);
                
        for g= 1: length(field_max_inds)-1
            for arena_count= g+1: length(field_max_inds)
                hyper_pt_distances(g,arena_count) = Distance(field_max_inds(g,1), field_max_inds(g,2),...
                    field_max_inds(arena_count,1), field_max_inds(arena_count,2));
            end
        end
        
        mean_distance_diff(count)= nanmean2(hyper_pt_distances);
        count=count+1;
    
    end
  
end

%end

mean_distance_diff