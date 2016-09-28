function MainRemappingData2
% Date: 17 of May 2015
dbstop if error

parms.dir_load_data = 'C:\Users\Dori\Desktop\saved_mat\saved_mat';
parms.dir_save_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\third';

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=5;
parms.sigma = 3;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;



% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    % find posx and posy
    
    for cell_num=1:length(cells)
        
        %fig=figure;
        
        mySpikes=[];
        
        zone_mats=[];
        
        gridness_scores=[];
        
        for arena_count= 1:length(tint);
            
            if tets(cell_num)==1
                mySpikes=find(cutTet1{arena_count}==cells(cell_num));
            elseif tets(cell_num)==2
                mySpikes=find(cutTet2{arena_count}==cells(cell_num));
            elseif tets(cell_num)==3
                mySpikes=find(cutTet3{arena_count}==cells(cell_num));
            elseif tets(cell_num)==4
                mySpikes=find(cutTet4{arena_count}==cells(cell_num));
            end
            
            if length(tint(arena_count).tetrode)>= tets(cell_num) %why is tetrode 4 missing when spikes are supposedly found at tetrode 4? file_name=216_6_6_12
                
                len= length(tint(arena_count).tetrode(tets(cell_num)).pos_sample); % why is mySpikes longer than the length of pos_sample???
                mySpikes(mySpikes>len)=[]; % why is mySpikes longer than the length of pos_sample???
                mySpikes=mySpikes';
                
                myPosSample=[];
                myPosSample= tint(arena_count).tetrode(tets(cell_num)).pos_sample(mySpikes);
                
                spk_x=[];
                spk_y=[];
                spk_x=tint(arena_count).pos.xy(myPosSample,1);
                spk_y=tint(arena_count).pos.xy(myPosSample,2);
                
                pos_t=[];
                pos_t=zeros([length(tint(arena_count).pos.xy), 1]);
                pos_t(1)=0.02;
                
                for h=2:length(pos_t);
                    pos_t(h)= pos_t(h-1) +0.02;
                end
                
                % create rate mat
                
                pos_x=[];
                pos_y=[];
                pos_x= tint(arena_count).pos.xy(:,1);
                pos_y= tint(arena_count).pos.xy(:,2);
                
                spk_t=[];
                
                for h=1:length(spk_x)
                    if isnan(spk_x(h))
                        spk_t(h)=NaN;
                    elseif ~isnan(spk_x(h))
                        spk_xy_ind= find(pos_x==spk_x(h) & pos_y==spk_y(h));
                        spk_t(h)= pos_t(spk_xy_ind(1));
                    end
                end
                
                [rate_mat]=CreateRateMap(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t,parms);
                
                S.spk_t=spk_t;
                S.pos_mean_x=pos_x;
                S.pos_mean_y=pos_y;
                S.pos_t=pos_t;
                S.spk_x= spk_x;
                S.spk_y=spk_y;
                
                S=get_zones(rate_mat,S);
                
                autocorr = Cross_Correlation(rate_mat,rate_mat);
                
                R_outer = FindROuter(autocorr,parms);
                [gridness2] = GridnessRadius(autocorr,parms,R_outer,i);
                
                gridness_scores(arena_count) = gridness2;
                
                if length(S.sorted_means) >=3;
                    
                    norm_third_max= ThirdPeakNormIndex(S.sorted_means, S.peak_rates, S.peak_zone_mat);
                    third_anchor_indices(arena_count,1)= norm_third_max(1);
                    third_anchor_indices(arena_count,2)= norm_third_max(2);
                    
                else
                    disp(sprintf('removed %s %d because not enough fields', file_name, cell_num));
                    third_anchor_indices(arena_count,1)= nan;
                    third_anchor_indices(arena_count,2)= nan;
                    break;
                end
                
                %  pivot_strengths{arena_count} = S.pivot;
                
                %else
                %  disp(sprintf('ERROR in %s', file_name));
                % end
                
                if arena_count == 1
                    rate_mats.arena1= rate_mat;
                elseif arena_count==2
                    rate_mats.arena2= rate_mat;
                elseif arena_count==3
                    rate_mats.arena3 = rate_mat;
                elseif arena_count==4
                    rate_mats.arena4= rate_mat;
                elseif arena_count==5
                    rate_mats.arena5= rate_mat;
                end
                
                
                if arena_count == 1
                    zone_mats.arena1= S.zone_mat;
                elseif arena_count==2
                    zone_mats.arena2= S.zone_mat;
                elseif arena_count==3
                    zone_mats.arena3= S.zone_mat;
                elseif arena_count==4
                    zone_mats.arena4= S.zone_mat;
                elseif arena_count==5
                    zone_mats.arena5= S.zone_mat;
                end
                
            end
            
            
        end
        %%THIRD
        
        if any(gridness_scores >= 0.3)
            
            third_anchor_pt_distances= nan(5);
            
            for g= 1: length(third_anchor_indices)-1
                for h= g+1: length(third_anchor_indices)
                    third_anchor_pt_distances(g,h) = Distance(third_anchor_indices(g,1), third_anchor_indices(g,2),...
                        third_anchor_indices(h,1), third_anchor_indices(h,2));
                end
            end
            
            
            
            
            
            distance_matrix = third_anchor_pt_distances;
            
            
            cluster_matrix= nan(5);
            
            for g= 1:length(distance_matrix)-1;
                for h= g+1:length(distance_matrix);
                    if distance_matrix(g,h) <= 0.35
                        cluster_matrix(g,h)= 1;
                    elseif distance_matrix(g,h) > 0.35
                        cluster_matrix(g,h)= 0;
                    end
                end
            end
            
            cluster_score= sum(sum(cluster_matrix==1))/ sum(sum(~isnan(cluster_matrix)));
            
            
            
            cd(parms.dir_save_data);
            save(sprintf('%s_cell%d.mat', file_name, cells(cell_num)), 'S', 'third_anchor_indices', 'third_anchor_pt_distances', 'rate_mats', 'zone_mats', 'cluster_matrix', 'cluster_score');
            cd(parms.dir_load_data);
            
        end
    end
    
    
    
end

close all;
clear anchor_indices;
clear second_anchor_indices;
clear min_anchor_indices;
clear third_anchor_indices;
clear S;

% clearvars -except i parms file_name file_names cell_num arena_count tint



end





