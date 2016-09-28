% Date: 20 of Jan 2016
dbstop if error

parms.dir_load_data = 'C:\Users\Dori\Desktop\saved_mat\saved_mat';

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=5;
parms.sigma = 3;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count=1;
count2=1;
norm_hyperinds_for_shuffling=[];
% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    % find posx and posy
    
    for cell_num=1:length(cells)
        
        mySpikes=[];
        
        gridness_scores= nan(1,length(tint));
        max_inds=cell(1,length(tint));
        arena_peak= cell(1,length(tint));
        
        rate_mats= cell(1,length(tint));
        zone_mats= cell(1,length(tint));
        max_indices= cell(1,length(tint));
        peak_rates_arena= cell(1,length(tint));
        
        hyper_norm_ranking= nan(1,length(tint));
        norm_ranks_all_arenas= cell(1,length(tint));
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
                
                myPosSample= tint(arena_count).tetrode(tets(cell_num)).pos_sample(mySpikes);
                
                spk_x=tint(arena_count).pos.xy(myPosSample,1);
                spk_y=tint(arena_count).pos.xy(myPosSample,2);
                
                pos_t=zeros([length(tint(arena_count).pos.xy), 1]);
                pos_t(1)=0.02;
                
                for h=2:length(pos_t);
                    pos_t(h)= pos_t(h-1) +0.02;
                end
                
                % create rate mat
                
                pos_x= tint(arena_count).pos.xy(:,1);
                pos_y= tint(arena_count).pos.xy(:,2);
                
                spk_t=nan(1,length(spk_x));
                
                for h=1:length(spk_x)
                    if isnan(spk_x(h))
                        spk_t(h)=NaN;
                    elseif ~isnan(spk_x(h))
                        spk_xy_ind= find(pos_x==spk_x(h) & pos_y==spk_y(h));
                        spk_t(h)= pos_t(spk_xy_ind(1));
                    end
                end
                
                [rate_mat]=CreateRateMap(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t,parms);
                
                autocorr = Cross_Correlation(rate_mat,rate_mat);
                
                R_outer = FindROuter(autocorr,parms);
                [gridness2] = GridnessRadius(autocorr,parms,R_outer,i);
                
                gridness_scores(arena_count) = gridness2;
                
                %save all info for future shuffling
                max_inds= FindMaxIndsRateMap(rate_mat);
                auto_max_inds= FindAutoMaxInds(autocorr);
                PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);
                max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, 1.8);
                
                arena_peak{arena_count}= size(rate_mat);
                
                rate_mats{arena_count}= rate_mat;
                max_indices{arena_count}= max_inds;
                
                
                peak_rates= nan(1,length(max_inds));
                for cen= 1:length(max_inds);
                    peak_rates(cen)= rate_mat(max_inds(cen,1), max_inds(cen,2));
                end
                               
                [zone_mats{arena_count},number_zone_mat]= CreateZoneMat(rate_mat, PF_radius, max_inds, peak_rates);
                
                peak_rates_arena{arena_count}= peak_rates;
                
                %Downscale spikes
                
                [spk_x_inds, spk_y_inds]= ConvertCoordinates(rate_mat, 5, spk_x, spk_y);
                [pos_x_inds, pos_y_inds]= ConvertCoordinates(rate_mat, 5, pos_x, pos_y);
                
                num_of_fields= length(max_inds);
                number_spks_per_field=zeros(1,num_of_fields);
                for h=1:length(spk_x);
                    if ~isnan(spk_x)
                        zone_num= number_zone_mat(spk_x_inds(h), spk_y_inds(h));
                        if zone_num ~= 0
                            number_spks_per_field(zone_num)= number_spks_per_field(zone_num)+1;
                        end
                    end
                end
                
                max_index=find(peak_rates==max(peak_rates));
                
                spikes_in_hyperfield=number_spks_per_field(max_index);
                
                [~, ~, ranking] = unique(number_spks_per_field);
                hyper_spike_rank= ranking(max_index);
                
                norm_ranking= ranking/ length(ranking);
                hyper_norm_ranking(arena_count)= norm_ranking(max_index);
                
                norm_ranks_all_arenas{arena_count}= norm_ranking;
                
            end
        end %end of arena_count
        
        hyper_inds_use= [];
        norm_hyper_inds= cell(1,length(max_inds));
        if any(gridness_scores >= 0.3) %&& any(hyper_norm_ranking < 0.5) %if gridness score >0.3 in at least 1 arena
            
            norm_ranks{count}=hyper_norm_ranking;
            
            count=count+1;
        end
        
        hyper_inds_use= [];
        norm_hyper_inds= cell(1,length(max_inds));
        if any(gridness_scores >= 0.3) && any(hyper_norm_ranking < 0.5) %if gridness score >0.3 in at least 1 arena
            
            norm_ranks_only_some{count2}=hyper_norm_ranking;
            
            count2=count2+1;
        end
        
        %             % save info if at least one arena of gridness >= 0.3
        %             all_rate_mats{count}= rate_mats;
        %             all_zone_mats{count}= zone_mats;
        %             all_max_inds{count}= max_indices;
        %             all_peak_rates{count}= peak_rates_arena;
        %
        %             for h= 1:length(max_indices)
        %                 hyper_inds_use= (max_indices{h});
        %                 size_use= arena_peak{h};
        %                 hyper_inds_use(:,1)=hyper_inds_use(:,1)/size_use(1);
        %                 hyper_inds_use(:,2)=hyper_inds_use(:,2)/size_use(2);
        %                 norm_hyper_inds{h}= hyper_inds_use;
        %             end
        %
        %             norm_hyperinds_for_shuffling{count}= norm_hyper_inds;
        %
        %             title=sprintf('%s cell %d', file_name, cell_num);
        %             filenames{count}= title;
        %
        %             count=count+1;
        %         end % end of if at least one 0.3 gridness score
        
    end %end of cell_num
    
end %end of file_names


save('norm ranks rescaling', 'norm_ranks', 'norm_ranks_only_some', 'hyperfield_norm_inds')

disp('');

%cd('N:\users\rebekkah\results and info of analysis');
%save('file names rescaling data', 'filenames');

% save('info for cluster score shuffling', 'norm_hyperinds_for_shuffling');
% save('rescaling arenas info', 'all_rates_mats', 'all_zone_mats', ...
%     'all_max_inds', 'all_peak_rates')