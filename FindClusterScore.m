% Shuffles cluster scores
load('rescaling arenas data info.mat','norm_hyperinds','peak_rates_all')

% shuffle_times=10000;
% 
% mean_cluster_scores_pt_45=nan(1,shuffle_times);
% mean_cluster_scores_pt_35=nan(1,shuffle_times);
% mean_cluster_scores_pt_25=nan(1,shuffle_times);
% 
% for k=1:shuffle_times
%     
%     cluster_scores_25= nan(1,length(norm_hyperinds));
%     cluster_scores_35= nan(1,length(norm_hyperinds));
%     cluster_scores_45= nan(1,length(norm_hyperinds));
    %remove above if not shuffling
    
    for i= 1:length(norm_hyperinds) %opens all cells
        
        peak_rates_all_arenas=peak_rates_all{i};
        hyper_inds=norm_hyperinds{i};
        
        len=length(peak_rates_all_arenas);
        
        hyperinds=nan(len,2);
        for arena_count=1:len    %opens all arenas
            peak_rates=peak_rates_all_arenas{arena_count};
            all_inds=hyper_inds{arena_count};
            
            %add for shuffling:
           % peak_rates=Shuffle(peak_rates);
            
            max_peak_ind= find(peak_rates==max(peak_rates));
            
            hyperinds(arena_count,:)= all_inds(max_peak_ind,:);
        end
        
        hyper_pt_distances= nan(len);
        % finds distances between all norm hyperfields in all arenas
        for g= 1: len-1
            for h= g+1:len
                hyper_pt_distances(g,h) = Distance(hyperinds(g,1), hyperinds(g,2),...
                    hyperinds(h,1), hyperinds(h,2));
            end
        end
        
        distance_matrix = hyper_pt_distances;
        
        cluster_matrix_25= nan(len);
        cluster_matrix_35= nan(len);
        cluster_matrix_45= nan(len);
        for g= 1:length(distance_matrix)-1;
            for h= g+1:length(distance_matrix);
                if distance_matrix(g,h) <= 0.45
                    cluster_matrix_45(g,h)= 1;
                elseif distance_matrix(g,h) > 0.45
                    cluster_matrix_45(g,h)= 0;
                end
                
                if distance_matrix(g,h) <= 0.35
                    cluster_matrix_35(g,h)= 1;
                elseif distance_matrix(g,h) > 0.35
                    cluster_matrix_35(g,h)= 0;
                end
                
                if distance_matrix(g,h) <= 0.25
                    cluster_matrix_25(g,h)= 1;
                elseif distance_matrix(g,h) > 0.25
                    cluster_matrix_25(g,h)= 0;
                end
            end
        end
        
        cluster_scores_25(i)= sum(sum(cluster_matrix_25==1))/ sum(sum(~isnan(cluster_matrix_25)));
        cluster_scores_35(i)= sum(sum(cluster_matrix_35==1))/ sum(sum(~isnan(cluster_matrix_35)));
        cluster_scores_45(i)= sum(sum(cluster_matrix_45==1))/ sum(sum(~isnan(cluster_matrix_45)));
        
    end

    save('cluster scores', 'cluster_scores_25','cluster_scores_35','cluster_scores_45')
 % remove below if not shuffling:   
%     mean_cluster_scores_pt_25(k)= mean(cluster_scores_25);
%     mean_cluster_scores_pt_35(k)= mean(cluster_scores_35);
%     mean_cluster_scores_pt_45(k)= mean(cluster_scores_45);
%     
%     k
%     
% end
% 
% save('cluster score shuffled results', 'mean_cluster_scores_pt_25', ...
%     'mean_cluster_scores_pt_35', 'mean_cluster_scores_pt_45')
