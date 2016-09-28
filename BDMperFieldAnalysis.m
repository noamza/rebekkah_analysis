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

BDMs_all_cells=cell(1,length(file_names)); %where answers will be stored for each cell

% enumerate on cells
for i= 1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    % if i~=5
    
    % calculate the the rat's position (using average of both leds)
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
    
    % dist_boundary=all_distances{i};
    
       spk_t=Cell.spk.t;
    
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
    
    BDMs=cell(1,length(max_inds)); %where all bursts will be
    for h= 1:length(spike_bursts)
        % BDMs_final=cell(1,length(spike_bursts{h}));
        SDM=nan(1,length(spike_bursts{h}));
        for b=1:length(spike_bursts{h});
            dist=nan(1,length(max_inds));
            for cen= 1:length(max_inds)
                spk_num_burst= spk_nums{h};
                spk_num=spk_num_burst(b);
                %finds closest field center to spike
                if ~isnan(spk_x(spk_num))
                    dist(cen)=Distance(spk_x(spk_num),spk_y(spk_num), max_inds(cen,1), max_inds(cen,2));
                end
            end
            
            if ~isnan(dist)
            min_ind=find(dist==min(dist));
            min_ind=min_ind(1); %in case two minimums
            SDM(b)=dist(min_ind); %finds the SDM of each spk in the spike burst
            end
        end
        
        SDM(isnan(SDM))=[];
        
        if length(SDM) > 1
        
        BDM= mean(SDM);
        BDM_final=BDMs{min_ind};
        BDM_final(end+1)=BDM/r; %divide by radius
        BDMs{min_ind}=BDM_final; 
        
        
        end
    end
    
    if length(SDM) >1
    BDMs_all_cells{i}=BDMs;
    
    hyper_ind=find(Cell.peak_rates==max(Cell.peak_rates));
    BDM_field(1)=mean(BDMs{hyper_ind});
    BDM_all(i,1)=BDM_field; 
    BDM_wo_hyper=BDMs;
    BDM_wo_hyper(hyper_ind)=[];
    BDM_all(i,2)= mean(cell2mat(BDM_wo_hyper));
    
    hyper_inds(i)= hyper_ind;
    end
end

disp('')

nanmean(BDM_all(:,1))
std(BDM_all(:,1))/sqrt(length((BDM_all(:,1))))
nanmean(BDM_all(:,2))
std(BDM_all(:,2))/sqrt(length((BDM_all(:,2))))

figure;
for h=1:length(BDMs_all_cells)
    if h~=52
        BDMs=BDMs_all_cells{h};
        BDMs{hyper_inds(h)}=[];
        for j=1:length(BDMs)
            len=length(BDMs);
            for k=1:length(BDMs{j})
                BDM=BDMs{j};
                plot((ones(1,length(BDM))*j)/len,BDM, 'o');
                hold on;
            end
        end
    end
end
 
for h=1:length(BDMs_all_cells)
     BDMs=BDMs_all_cells{h};
    for k=1:length(BDMs{hyper_inds})
            BDM=BDMs{hyper_inds(h)}
            plot((ones(1,length(BDM))*j)/len,BDM, 'ro'); 
            hold on;
    end
end
        
    
   

disp('')
