cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat', 'zone_mats_ns_all',...
    'norm_max_index', 'max_indices_ns', 'rate_mats_all', 'peak_rates_sm_all','max_indices_sm')
load('rat date info.mat')

 means=nan(1,1000);
%
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
                
                
                strfind1=(strfind(rats_use,all_rats{c}));
                strfind2=strfind(dates_use,all_dates{b});
                strfind3=strfind(sesh_use,all_seshs{a});
                
                inds1=[];
                for h=1:length(strfind1)
                    if strfind1{h}==1 & strfind2{h}==1 & strfind3{h}==1
                        inds1(end+1)= h;
                    end
                end
                
                n=3; m=3;
                if length(inds1)>=3
                    
                    %                        fig=figure;
                    %                        for h=1:length(inds1)
                    %
                    %                        subplot(n,m,h)
                    %                        imagesc(zone_mats_ns_all{inds1(h)})
                    %                        axis off; axis image;
                    %
                    %                        end
                    %
                    %
                    %                      %   for shuffling: 
                                             rands=nan(1,length(inds1));
                                             for h=1:length(inds1)
                                                 ran= Shuffle(1:length(max_indices_ns{inds1(h)}));
                    rands(h)=ran(1);
                                             end
                    
                    
                    hyper_pt_distances=nan(length(inds1),length(inds1));
                    for g=1:length(inds1)-1
                        for h= g+1:length(inds1)
                            ind3=inds1(g);
                            ind2=inds1(h);
                            
                            max_indices1= max_indices_ns{ind3};
                            max_indices2= max_indices_ns{ind2};
                            
                            
                            %                                % for shuffling:
                                                             hyper_ind1=max_indices1(rands(g),:)./size(zone_mats_ns_all{ind3});
                                                            hyper_ind2=max_indices2(rands(h),:)./size(zone_mats_ns_all{ind2});
                            
%                             hyper_ind1=norm_max_index(ind3,:);
%                             hyper_ind2=norm_max_index(ind2,:);
                            
                            
                            hyper_pt_distances(g,h) = Distance(hyper_ind1(1), hyper_ind1(2),...
                                hyper_ind2(1), hyper_ind2(2));
                        end
                    end
                    
                    cluster_score= sum(hyper_pt_distances(:)<=0.35)/sum(~isnan(hyper_pt_distances(:)));
                    
                    clusters(count)=cluster_score;
                    
                    
%                                             subplot(n,m,length(inds1))
%                                             title(sprintf('%0.2f', cluster_score));
%                     %
%                     %
%                         cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7\simutaneous images')
%                                             saveas(fig, sprintf('%d.jpg', count));
%                                             count=count+1;
                    
                  
                    
%                     
%                     if count == 7 | count==10 | count==36
%                         
%                         letter={'B', 'A', 'C'};
%                         
%                         if count==10
%                             step=[0 0 0 0 0];
%                             com=6;
%                         elseif count==7
%                             step=[6 6 6];
%                             com=12;
%                         elseif count==36
%                             step=[12 12 12 12 12 13 13];
%                             com=24;
%                         end
%                         
%                       
%                        % PrintSimutaneousImage(max_indices_sm,letter, inds1, rate_mats_all, step, com, rats_use,dates_use, sesh_use, peak_rates_sm_all)
%                           
%                         
%                         disp('')
%                     end
                    
                    % close all
                    count=count+1;
                end
                
                
                
                disp('')
                
                
            end
        end
    end
end

% subplot(6,6, [25 26 27])
% hist(clusters);
% xlabel('Cluster score')
% ylabel('# of cells')
% text(-0.25,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% 
% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 23.3], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 23.3])


disp('')
  means(k)=mean(clusters);

save('simutaniety shuffled data', 'means')
     k
 end
%
%
 figure; hist(means)