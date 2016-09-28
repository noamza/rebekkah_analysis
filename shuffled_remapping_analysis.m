parms.dir_load_data = 'C:\Users\Dori\Desktop\saved_mat\saved_mat';

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=5;
parms.sigma = 3;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};



for shuffle=1:250;
    
    count=0;
    
    % enumerate on cells
    for i =1:length(file_names)
        cd(parms.dir_load_data);
        file_name = file_names{i};
        load(file_name);
        
        % find posx and posy
        
        
        
        
        
        for cell_num=1:length(cells)
            
            
            S.rate_mat=[];
            
            S.zone_mat=[];
            
            mySpikes=[];
            
            pivot_strengths=[];
            
            gridness2=[];
            
            gridness=[];
            
            anchor_indices=[];
            
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
                    
                    
                    
                    %shuffles results
                    peak_rates_list = S.peak_rates(randperm(length(S.peak_rates)));
                    
                    zone_num = find(peak_rates_list==S.sorted_means(end));
                    
                    [row, col] = find(S.peak_zone_mat == zone_num);
                    
                    index = [];
                    index(1) = row;
                    index(2) = col;
                    
                    [size_x, size_y] = size(S.peak_zone_mat);
                    
                    norm_max_index(1)= index(1)/size_x;
                    norm_max_index(2)=index(2)/size_y;
                    % results are shuffled above.
                    
                    
                    
                    % indices of shuffled max peaks
                    anchor_indices(arena_count, 1)= norm_max_index(1);
                    anchor_indices(arena_count, 2)= norm_max_index(2);
                    
                    
                    
                    
                    
                    % finds distances between all points
                    
                    
                    
                    %         if any(gridness_scores >= 0.3)
                    %
                    %             anchor_pt_distances= nan(5);
                    %
                    %             for g= 1: length(anchor_indices)-1
                    %                 for h= g+1: length(anchor_indices)
                    %                     anchor_pt_distances(g,h) = Distance(anchor_indices(g,1), anchor_indices(g,2),...
                    %                         anchor_indices(h,1), anchor_indices(h,2));
                    %                 end
                    %             end
                    %
                    %         end
                    
                    
                    
                    
                    
                    
                    
                    close all;
                    
                    clear S;
                    
                    
                    
                    % clearvars -except i parms file_name file_names cell_num arena_count tint
                    
                    
                    
                    
                    
                end
                
                
            end
            
            
            
            
            
        
        
        
        above_grid_scores= (gridness_scores>= 0.3);
        
        if sum(above_grid_scores) >= 1
            
            anchor_pt_distances= nan(5);
            
            for g= 1: length(anchor_indices)-1
                for h= g+1: length(anchor_indices)
                    anchor_pt_distances(g,h) = Distance(anchor_indices(g,1), anchor_indices(g,2),...
                        anchor_indices(h,1), anchor_indices(h,2));
                end
            end
            
            
            distance_matrix = anchor_pt_distances;
            
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
            
            count=count+1;
            cluster_scores(count) = cluster_score;
            
        end       
            
        end
        
        
        
    end
    
    num_cluster=length(cluster_scores);
    
    num_above_three(shuffle)= sum(cluster_scores >= 0.3); %number of cells with cluster scores above 0.3
    num_above_four(shuffle)=sum(cluster_scores >= 0.4); %number of cells with cluster scores above 0.4
    num_above_six(shuffle)=sum(cluster_scores >= 0.6); %number of cells with cluster scores above 0.6
    mean_cluster(shuffle)=mean(cluster_scores); %the mean cluster score of each shuffled set
   
    
    
    
end

save('shuffed results FULL DATA 2', 'num_above_three', 'num_above_four', 'num_above_six', 'mean_cluster')
