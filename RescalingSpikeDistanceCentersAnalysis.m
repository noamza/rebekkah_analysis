function RescalingSpikeDistanceCentersAnalysis
% Date: 17 of May 2015
dbstop if error

parms.dir_load_data2= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';
parms.dir_load_data = 'C:\Users\Dori\Desktop\saved_mat\saved_mat';
parms.dir_save_data2 = 'N:\users\rebekkah\are spks closer to max peaks\barry jeffery results';


parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=5;
parms.sigma = 3;

dir_name=parms.dir_load_data2;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;



% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    file_name2=file_name(1:end-10);
    load(file_name2);
    
    % find posx and posy
    
   
    for cell_num=1:length(cells)
        
        
        
        mySpikes=[];
       
       
        
        for arena_count= 1:length(tint);
            
             ordered_distances=[];
            
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
                
                [ordered_distances,  peak_score]= SpikeDistanceFromCenter(S.max_inds, spk_x, spk_y, S.peak_rates, S.sorted_means, pos_x, pos_y);

                max_field_spike_distance(arena_count) = ordered_distances(end);
                non_max_spike_mean_distance(arena_count)= nanmean(ordered_distances(1:end-1));
                
                
          
      
        
        
        
    end
    
    
    
        end

  cd(parms.dir_save_data2);
            save(sprintf('%s_cell%d.mat', file_name, cells(cell_num)), 'max_field_spike_distance', 'non_max_spike_mean_distance', 'peak_score');
            cd(parms.dir_load_data);
            
            figure; scatter(max_field_spike_distance, non_max_spike_mean_distance); hold on;
            title(sprintf('%0.2f',peak_score)); 

x= 0:30
y=0:30

plot(x,y,'-')


disp('')



    end
end

            