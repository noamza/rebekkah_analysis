% Shuffles cluster scores

load('rescaling arenas data info.mat', 'all_norm_inds', 'peak_rates_all')

 mean_cluster_scores_pt_45=nan(1,1000);
 mean_cluster_scores_pt_35=nan(1,1000);
 mean_cluster_scores_pt_25=nan(1,1000);
% 
 for k=1:10000
%     
    cluster_scores_25= nan(1,length(all_norm_inds));
    cluster_scores_35= nan(1,length(all_norm_inds));
    cluster_scores_45= nan(1,length(all_norm_inds));
    
    for i= 1:length(all_norm_inds) %opens all cells
        
        hyper_pt_distances= nan(5);
        max_inds= all_norm_inds{i};
        peak_rates_all_arenas= peak_rates_all{i};
        
        hyper_max_inds=nan(length(max_inds),2);   
        for arena_count=1:length(max_inds)    %opens all arenas
            
            peak_rates=peak_rates_all_arenas{arena_count};
            
            %peak_rates(peak_rates==max(peak_rates))=nan;
            max_peak_ind= find(peak_rates==max(peak_rates)); %changed max to min
            
            shuffle_inds= randperm(length(max_inds{arena_count}));
            inds=max_inds{arena_count};     %all the norm inds in arena 
          %   shuffles data
             inds= inds(shuffle_inds,:);
            
            hyper_max_inds(arena_count,:)= inds(max_peak_ind,:);
        end
        
        for g= 1: length(hyper_max_inds)-1
            for arena_count= g+1: length(hyper_max_inds)
                hyper_pt_distances(g,arena_count) = Distance(hyper_max_inds(g,1), hyper_max_inds(g,2),...
                    hyper_max_inds(arena_count,1), hyper_max_inds(arena_count,2));
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
        
        cluster_matrix_25= nan(5);
        
        for g= 1:length(distance_matrix)-1;
            for arena_count= g+1:length(distance_matrix);
                if distance_matrix(g,arena_count) <= 0.25
                    cluster_matrix_25(g,arena_count)= 1;
                elseif distance_matrix(g,arena_count) > 0.25
                    cluster_matrix_25(g,arena_count)= 0;
                end
            end
        end
        
        cluster_matrix_45= nan(5);
        
        for g= 1:length(distance_matrix)-1;
            for arena_count= g+1:length(distance_matrix);
                if distance_matrix(g,arena_count) <= 0.45
                    cluster_matrix_45(g,arena_count)= 1;
                elseif distance_matrix(g,arena_count) > 0.45
                    cluster_matrix_45(g,arena_count)= 0;
                end
            end
        end
        
        cluster_scores_25(i)= sum(sum(cluster_matrix_25==1))/ sum(sum(~isnan(cluster_matrix_25)));
        cluster_scores_35(i)= sum(sum(cluster_matrix==1))/ sum(sum(~isnan(cluster_matrix)));
        cluster_scores_45(i)= sum(sum(cluster_matrix_45==1))/ sum(sum(~isnan(cluster_matrix_45)));
        
    end
    
    % saves data set Cluster score (uncomment shuffling info for this)
    % save('Cluster scores second max', 'cluster_scores_25', 'cluster_scores_35', 'cluster_scores_45');
    
     mean_cluster_scores_pt_25(k)= mean(cluster_scores_25);
    mean_cluster_scores_pt_35(k)= mean(cluster_scores_35);
     mean_cluster_scores_pt_45(k)= mean(cluster_scores_45);
%     
     k
      end
% 
 save('cluster score shuffled results', 'mean_cluster_scores_pt_25', ...
     'mean_cluster_scores_pt_35', 'mean_cluster_scores_pt_45')
