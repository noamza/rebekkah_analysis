dbstop if error

% Shuffles cluster scores

load('rescaling arenas info.mat', 'all_norm_inds', 'all_peak_rates')

mean_cluster_scores=nan(1,5000);
median_cluster_scores=nan(1,5000);

for k=1:5000
    %
    cluster_scores_25= nan(1,length(all_norm_inds));
    
    for i= 1:length(all_norm_inds) %opens all cells
        
        hyper_pt_distances= nan(5);
        max_inds_all_arenas= all_norm_inds{i};
        peak_rates_all_arenas= all_peak_rates{i};
        
        % hyper_max_inds=nan(length(max_inds),2);
        rand_num= randperm(length(max_inds_all_arenas{1}));
        %rand_num= rand_num(1);
        
        peak_rates_orig= peak_rates_all_arenas{1};
        peak= peak_rates_orig(rand_num(1));
        
        max_peak= max(peak_rates_orig);
        %         peak=max_peak;
        
        %cannot be max peak:
        if peak == max_peak
            peak= peak_rates_orig(rand_num(2));
        end
        
        field_max_inds=nan(length(max_inds_all_arenas),2);
        for arena_count=2:length(max_inds_all_arenas)    %opens all arenas
            
            peak_rates=peak_rates_all_arenas{arena_count};
            
            %find ind of closest firing rate
            peak_rates_minus_peak= abs(peak_rates-peak);
            ind= find(peak_rates_minus_peak==min(peak_rates_minus_peak));
            
            max_inds= max_inds_all_arenas{arena_count};
            
            field_max_inds(arena_count,:)= max_inds(ind,:);
            
        end
        
        for g= 1: length(field_max_inds)-1
            for arena_count= g+1: length(field_max_inds)
                hyper_pt_distances(g,arena_count) = Distance(field_max_inds(g,1), field_max_inds(g,2),...
                    field_max_inds(arena_count,1), field_max_inds(arena_count,2));
            end
        end
        
        distance_matrix = hyper_pt_distances;
        cluster_matrix= nan(5);
        
        for g= 1:length(distance_matrix)-1;
            for arena_count= g+1:length(distance_matrix);
                if distance_matrix(g,arena_count) <= 0.35
                    cluster_matrix(g,arena_count)= 1;
                elseif distance_matrix(g,arena_count) > 0.35
                    cluster_matrix(g,arena_count)= 0;
                end
            end
        end
        
        cluster_scores(i)= sum(sum(cluster_matrix==1))/ sum(sum(~isnan(cluster_matrix)));
        
        
    end 
    
    % saves data set Cluster score (uncomment shuffling info for this)
    %     % save('Cluster scores updated', 'cluster_scores_25', 'cluster_scores_35', 'cluster_scores_45');
    %
    mean_cluster_scores(k)= mean(cluster_scores);
    
    k
    %
end
%
save('firing field cluster score shuffled results cant be hyperfield', 'mean_cluster_scores' ...
    )
