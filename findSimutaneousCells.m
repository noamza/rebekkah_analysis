cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat', 'zone_mats_ns_all', 'norm_max_index', 'max_indices_ns')
load('rat date info.mat')

means=nan(1,1000);

for k=1:1000
    
    count=1;
    for d=1:length(unique(study));
        
        inds=find(study==d);
        rats_use=rat(inds);
        dates_use=date(inds);
        sesh_use=session(inds);
        
        for c=1:length(unique(rats_use));
            
            all_rats=unique(rats_use);
            
            for b=1:length(unique(dates_use));
                
                all_dates=unique(dates_use);
                
                for a=1:length(unique(sesh_use));
                    
                    all_seshs=unique(sesh_use);
                    
                    
                    strfind1=strcmp(rats_use,all_rats{c});
                    strfind2=strcmp(dates_use,all_dates{b});
                    strfind3=strcmp(sesh_use,all_seshs{a});
                    
                    all=strfind1+strfind2+strfind3;
                    
                    sim_cell_inds=find(all==3);
                    
                    n=3; m=3;
                    if length(sim_cell_inds)>=3
                        %for individual images add:
%                                                fig=figure;
%                                                for h=1:length(sim_cell_inds)
%                         
%                                                subplot(n,m,h)
%                                                imagesc(zone_mats_ns_all{sim_cell_inds(h)})
%                                                axis off; axis image;
%                         
%                                                end
                        %
                        
                        %for shuffling:
                        rands=nan(1,length(sim_cell_inds));
                        for h=1:length(sim_cell_inds)
                            ran= Shuffle(1:length(max_indices_ns{sim_cell_inds(h)}));
                            rands(h)=ran(1);
                        end
                        
                        
                        hyper_pt_distances=nan(length(sim_cell_inds),length(sim_cell_inds));
                        for g=1:length(sim_cell_inds)-1
                            for h= g+1:length(sim_cell_inds)
                                ind3=sim_cell_inds(g);
                                ind2=sim_cell_inds(h);
                                
                                max_indices1= max_indices_ns{ind3};
                                max_indices2= max_indices_ns{ind2};
                                
                                
                                % for shuffling:
                                hyper_ind1=max_indices1(rands(g),:)./size(zone_mats_ns_all{ind3});
                                hyper_ind2=max_indices2(rands(h),:)./size(zone_mats_ns_all{ind2});
%                                 
                                % for non-shuffling:
%                                 hyper_ind1= norm_max_index(ind3,:);
%                                 hyper_ind2= norm_max_index(ind2,:);
                                
                                hyper_pt_distances(g,h) = Distance(hyper_ind1(1), hyper_ind1(2),...
                                    hyper_ind2(1), hyper_ind2(2));
                            end
                        end
                        
                        cluster_score= sum(hyper_pt_distances(:)<=0.35)/sum(~isnan(hyper_pt_distances(:)));
                        
                        clusters(count)=cluster_score;
                        
                          %for individual images add:
%                                                subplot(n,m,length(sim_cell_inds))
%                                                title(sprintf('%0.1f', cluster_score));
%                         
%                         
%                                                cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7\simutaneous images')
%                                                saveas(fig, sprintf('%d.jpg', count));
                                               
                        
                        count=count+1;
                        close all
                    end
                    
                    
                    
                    disp('')
                    
                  
                    
                end
            end
        end
    end
    
    
%       mean(clusters)
%       std(clusters)/sqrt(length(clusters)) 
      
      
    means(k)=mean(clusters);
    
    save('simutaniety shuffled data','means')
    k
end


figure; hist(means)