parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
%parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\data sets\images for adaptation analysis nonsmooth';
%parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\data sets\results for adaptation analysis nonsmooth';

parms.bin_size=3;
parms.sigma = 1.5;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

all_count=1;
% enumerate on cells
for i= 1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    if i~=5
        
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
        
        %finds the SDM for each spk
        %     SDM= nan(1,length(spk_x));
        %     for h=1:length(spk_x)
        %
        %         dist=nan(1,length(max_inds));
        %         if ~isnan(spk_x(h))    %if spk isn't nan
        %             for cen= 1:length(max_inds)
        %                 dist(cen)=Distance(spk_x(h),spk_y(h), max_inds(cen,1), max_inds(cen,2));
        %             end
        %             min_ind=find(dist==min(dist));
        %             min_ind=min_ind(1); %in case two minimums
        %             SDM(h)=Distance(spk_x(h),spk_y(h), max_inds(min_ind,1), max_inds(min_ind,2));
        %         end
        %     end
        %
        r= findPlaceFieldRadius(autocorr, auto_max_inds);
        %     SDM=SDM/r; %divide by field radius so can be compared to all cells
        
        %finds distance from border at each timepoint to define in boundary and
        %out of boundary for trajectory detection
        size_x= [min(pos_mean_x) max(pos_mean_x)];
        size_y= [min(pos_mean_y) max(pos_mean_y)];
        dist_boundary=nan(1,length(Cell.pos.t));
        for h=1:length(Cell.pos.t)
            [~,dist_boundary(h)]= findDistPtToBorder(size_x, size_y, [pos_mean_x(h) pos_mean_y(h)]);
        end
        
        % doesn't work, find finds index off by one
        %     inds= 2:length(dist_boundary)-1;
        %     ex= find(dist_boundary(inds) >= 0.1 & (dist_boundary(inds-1) <0.1 ));
        %     en= find(dist_boundary(inds) < 0.1 & (dist_boundary(inds-1) >=0.1 ));
        
        %finds indices of point of border exit and border entrance
        %(i.e. trajectory start and end)
        ex=[];
        en=[];
        count=1;
        count2=1;
        for inds= 2:length(dist_boundary)-1;
            if (dist_boundary(inds) >= 0.1 & (dist_boundary(inds-1) <0.1 ));
                ex(count)= inds;
                count=count+1;
            elseif (dist_boundary(inds) < 0.1 & (dist_boundary(inds-1) >=0.1 ));
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
            
            time_since_boundary= 0:0.02:length(traject_pos_t{traject})*0.02;
            time_since_boundary(end)= [];
            
            pos_t= traject_pos_t{traject};
            
            %find the spike positions where spks occur
            
            time_start= find(Cell.spk.t >= min(pos_t));
            
            if length(time_start) >1
                
                
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
                
                if length(timestamps) >=2  %if at least 2 spikes in trajectory
                    
                    traj_spk_x=spk_x(time_start:time_end);
                    traj_spk_y=spk_y(time_start:time_end);
                    
                    traj_len=nan(1,length(traj_spk_x));
                    SDM_traject=nan(1,length(traj_spk_x));
                    for h=1:length(traj_spk_x);
                        dist=nan(1,length(max_inds));
                        for cen= 1:length(max_inds)
                            dist(cen)=Distance(traj_spk_x(h),traj_spk_y(h), max_inds(cen,1), max_inds(cen,2));
                        end
                        min_ind=find(dist==min(dist));
                        min_ind=min_ind(1); %in case two minimums
                        SDM_traject(h)=Distance(traj_spk_x(h),traj_spk_y(h), max_inds(min_ind,1), max_inds(min_ind,2));
                        traj_len(h)= time_since_boundary(timestamps(h));
                    end
                    
                    SDM_traject=SDM_traject/r; %divide by field radius
                    
                    %remove nans and places where are nans
                    %                 nan_inds= find(isnan(SDM_traject));
                    %                 time_since_boundary(nan_inds)=[];
                    %                 SDM_traject(nan_inds)=[];
                    
                    %figure; scatter(time_since_boundary,SDM_traject); lsline;
                    slope=polyfit(traj_len,SDM_traject,1);
                    
                    slopes(all_count)=slope(1);
                    all_count=all_count+1;
                    
                end
            end
            end
            
        end %end of for all_traject
        
        i
        
    end %if not 5th file (corrupted)
end %end of for all files

pos_slopes= slopes;
pos_slopes(pos_slopes>0)=1;
pos_slopes(pos_slopes<=0)=0;

figure; hist(slopes);
mean(slopes)

for h=1:0.05:max(traj_lens)
    percent_pos(h/0.05)= sum(pos_slopes(traj_lens>=h)==1)/ length(pos_slopes(traj_lens>=h));
    x(h/0.05)=h;
end

figure; plot(x, percent_pos)

disp('')
