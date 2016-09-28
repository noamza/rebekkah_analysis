dbstop if error

% opens each cell file and cuts the entire session trajectory where the rat 
% exits and re-enters a boundary area (defined as 11 cm here)
% finds spike bursts for each trajectory (at least 3 bursts and at least
% 2 secs long) 
% finds slope of the mean distance of the bursts from the field center as a
% function of the time since it exited the boundary area

cd('N:\users\rebekkah\results and info of analysis')
load('err accumulation info at least 2 secs long at least 3 spike bursts.mat', 'all_distances')

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
cd(parms.dir_load_data)

parms.bin_size=3;
parms.sigma = 1.5;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

%all_distances=cell(1,length(file_names));

all_count=1;
% enumerate on cells
for i= 1:16 %length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    % if i~=5
    
    % calculate the the rat's head direction (using average of both leds)
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    %Giocomo defined radius differently for each field, will stick to our
    %definition
    autocorr=Cross_Correlation(Cell.rate_mat,Cell.rate_mat);
    auto_max_inds = FindAutoMaxInds(autocorr);
    
    %convert x and y coordinated to matrix indices used in rate map
    [pos_i, pos_j]= ConvertCoordinates(Cell.rate_mat, parms.bin_size,pos_mean_x,pos_mean_y);
    [spk_i, spk_j]= ConvertCoordinates(Cell.rate_mat, parms.bin_size,spk_x,spk_y);
    
    %convert matrix indices to x y coordinates
    inds_x_y= nan(length(Cell.max_inds),2);
    inds_x_y(:,1)= Cell.max_inds(:,2)* parms.bin_size + min(pos_mean_x);
    inds_x_y(:,2)= Cell.max_inds(:,1)* parms.bin_size + min(pos_mean_y);
    
    max_inds=inds_x_y;
     
    r= findPlaceFieldRadius(autocorr, auto_max_inds);
    
    %finds distance from border at each timepoint to define in boundary and
    %out of boundary for trajectory detection
%     size_x= [min(pos_mean_x) max(pos_mean_x)];
%     size_y= [min(pos_mean_y) max(pos_mean_y)];
%     dist_boundary=nan(1,length(Cell.pos.t));
%     for h=1:length(Cell.pos.t)
%         [dist_boundary(h),~]= findDistPtToBorder(size_x, size_y, [pos_mean_x(h) pos_mean_y(h)]);
%     end
%     
%     all_distances{i}=dist_boundary;
    
dist_boundary=all_distances{i}; 

    % doesn't work, find finds index off by one
    %     inds= 2:length(dist_boundary)-1;
    %     ex= find(dist_boundary(inds) >= 0.1 & (dist_boundary(inds-1) <0.1 ));
    %     en= find(dist_boundary(inds) < 0.1 & (dist_boundary(inds-1) >=0.1 ));
    
    %finds indices of point of border exit and border entrance
    %(i.e. trajectory start and end)
    %use 11 cm for boudary threshold (same as Giocomo paper)
    ex=[];
    en=[];
    count=1;
    count2=1;
    for inds= 2:length(dist_boundary)-1;
        if (dist_boundary(inds) >= 11 & (dist_boundary(inds-1) < 11 ));
            ex(count)= inds;
            count=count+1;
        elseif (dist_boundary(inds) < 11 & (dist_boundary(inds-1) >= 11 ));
            en(count2)= inds;
            count2=count2+1;
        end
    end
    
    en(en < ex(1))=[]; %disregard if enters boundary before exiting
    
    % cuts trajectories and saves pos_t
    lenn=min(length(ex), length(en));
    traject_pos_t=cell(1,length(lenn));
    traject_pos_x=cell(1,length(lenn));
    traject_pos_y=cell(1,length(lenn));
    for h=1:lenn
        traject_pos_t{h}= Cell.pos.t(ex(h):en(h)-1);
        traject_pos_x{h}= pos_mean_x(ex(h):en(h)-1);
        traject_pos_y{h}= pos_mean_y(ex(h):en(h)-1);
    end
    
    for traject=1:length(traject_pos_t);
        
        if length(traject_pos_t{traject}) > 100 % only look at trajectories > 2sec
            
            time_since_boundary= 0:0.02:length(traject_pos_t{traject})*0.02;
            time_since_boundary(end)= [];
            
            pos_t= traject_pos_t{traject};
            
            %find the spike positions where spks occur
            
            time_start= find(Cell.spk.t >= min(pos_t));
            
            if length(time_start) >= 1
                
                time_start=time_start(1);
                time_end= find(Cell.spk.t <= max(pos_t));
                
                if length(time_end) >= 1
                    
                    time_end= time_end(end);
                    len= length(Cell.spk.t(time_start:time_end));
                    timestamps=nan(1,len);
                    for h=1:len
                        timestamp= find(pos_t >= Cell.spk.t(time_start+h-1));
                        timestamps(h)=timestamp(1);
                    end
                    
                    if length(timestamps) >=4  %if at least 4 spikes in trajectory (minimum nesseary for at least 2 spike bursts)
                        
                        spk_t= Cell.spk.t(time_start:time_end);
                        traj_spk_x=spk_x(time_start:time_end);
                        traj_spk_y=spk_y(time_start:time_end);
                        
                        count=1;
                        burst_count=1;
                        burst=[];
                        spike_bursts=[];
                        spk_num=[];
                        spk_nums=[];
                        for h=1:length(spk_t)-1
                            time_btwn= spk_t(h+1)-spk_t(h);
                            if time_btwn<=0.02
                                burst(burst_count:burst_count+1)=spk_t(h:h+1);
                                spk_num(burst_count:burst_count+1)=h:h+1;
                                burst_count=burst_count+1;
                            elseif time_btwn >0.02 & length(burst)> 1
                                spike_bursts{count}=burst;
                                spk_nums{count}=spk_num;
                                count=count+1;
                                burst=[];
                                spk_num=[];
                                burst_count=1;
                            end
                            
                            if h==length(spk_t)-1 & length(burst)>1
                                spike_bursts{count}=burst;
                                spk_nums{count}=spk_num;
                            end
                        end
                        
                        
                        if length(spike_bursts)>=2 %at least 3 spike bursts in trajectory
                            
                            BDM_traject=nan(1,length(spike_bursts));
                            traj_len=nan(1,length(spike_bursts));
                            for b=1:length(spike_bursts)
                                
                                spk_len=nan(1,length(spike_bursts{b}));
                                SDM=nan(1,length(spike_bursts{b}));
                                for h=1:length(spike_bursts{b});
                                    dist=nan(1,length(max_inds));
                                    for cen= 1:length(max_inds)
                                        x=spk_nums{b};
                                        dist(cen)=Distance(traj_spk_x(x(h)),traj_spk_y(x(h)), max_inds(cen,1), max_inds(cen,2));
                                    end
                                    min_ind=find(dist==min(dist));
                                    min_ind=min_ind(1); %in case two minimums
                                    SDM(h)=dist(min_ind); %finds the SDM of each spk in the spike burst
                                    
                                end
                                
                                
                                time=timestamps(spk_nums{b});
                                len=(time_since_boundary(time));
                                traj_len(b)=mean(len);
                                BDM_traject(b)= mean(SDM);
                                
                                
                            end
                            
                            BDM_traject=BDM_traject/r ; %divide by field radius
                            
                            if range(traj_len) ~= 0
                            
                            %remove nans and places where are nans
                            %                 nan_inds= find(isnan(SDM_traject));
                            %                 time_since_boundary(nan_inds)=[];
                            %                 SDM_traject(nan_inds)=[];
                            
                            %figure; scatter(time_since_boundary,SDM_traject); lsline;
                            slope=polyfit(traj_len,BDM_traject,1);
                            
                            slopes(all_count)=slope(1);
                            len_traj(all_count)= length(traject_pos_t{traject});                            
                            file(all_count)= i;
                            trajectory(all_count)=traject;
                            time_of_traject{all_count}=traj_len;
                            BDMs{all_count}=BDM_traject;
                            arena_sizes{all_count}=size(Cell.rate_mat);
                            all_count=all_count+1;
                            
%                             if i==any(1:16) 
%                                directory(all_count)='S';
%                             elseif i==any(17:21)
%                                 directory(all_count)='B';
%                             else
%                                 directory(all_count)='D';
%                             end
                            
                            
                            end 
                        end %end of spike bursts
                   end
                end
            end
            
        end %end if traject > 1 sec long 
    end %end of for all_traject
    
    i
    
    %end %if not 5th file (corrupted)
end %end of for all files

save('err accumulation info at least 2 secs long at least 2 spike bursts', ...
    'slopes', 'len_traj', 'file', 'all_distances', 'trajectory', 'time_of_traject',...
    'BDMs', 'arena_sizes')

pos_slopes= slopes;
pos_slopes(pos_slopes>0)=1;
pos_slopes(pos_slopes<=0)=0;

figure; hist(slopes);
mean(slopes)

count=1;
for h=1:0.05:max(len_traj)
    percent_pos(count)= sum(pos_slopes(len_traj>=h)==1)/ length(pos_slopes(len_traj>=h));
    x(count)=h;
    count=count+1;
end

figure; plot(x*0.02, percent_pos)

figure;
for h=1:length(time_of_traject)
   plot(time_of_traject{h},BDMs{h},'o');
   hold on;
   
end




disp('')
