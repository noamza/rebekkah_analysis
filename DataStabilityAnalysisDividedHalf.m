dbstop if error

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\3 datasets G3MD25';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

%info for rate map creation
parms.bin_size= 3;
parms.sigma= 1.5;

files_len= 85;

peak_rates_b=cell(1,files_len);
peak_rates_e=cell(1,files_len);

for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    if i ~= 5 %fifth file corrupted
        
        begin_inds=[];
        end_inds=[];
        
        pos_t_len= length(Cell.pos.t);
        
        %divides pos_t in half
        if mod(pos_t_len,2) ==0   % even number
            begin_inds= 1:pos_t_len/2;
            end_inds= (pos_t_len/2)+1:pos_t_len;
        elseif mod(pos_t_len,2) == 1    % odd number
            begin_inds= 1:(pos_t_len/2)+0.5;
            end_inds= (pos_t_len/2)+1.5:pos_t_len;
        end
        
        pos_t_begin= Cell.pos.t(begin_inds);   %first half of session
        pos_t_end= Cell.pos.t(end_inds); %second half of session
        
        tmp_begin = abs(Cell.spk.t-pos_t_begin(end)); %subtract last timestamp of first session
        tmp_end= abs(Cell.spk.t- pos_t_end(1)) ;% subtract first timestamp of second session
        [~, idx_begin] = min(tmp_begin); %index of first session
        [~, idx_end] = min(tmp_end) ;%index of first session
        spk_t_begin = Cell.spk.t(1:idx_begin); %spike vector of first sesh
        spk_t_end = Cell.spk.t(idx_end+1:length(Cell.spk.t)); %spike vector of second sesh
        
        % calculate the the rat's position (using average of both leds)
        pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
        % divide into half
        pos_mean_x_begin= pos_mean_x(begin_inds); %x position for first half of session
        pos_mean_x_end= pos_mean_x(end_inds); %x position for second half of session
        pos_mean_y_begin= pos_mean_y(begin_inds); %y position for first half of session
        pos_mean_y_end= pos_mean_y(end_inds); %y position for second half of session
        
        spk_x_begin= spk_x(1:idx_begin);
        spk_x_end= spk_x(idx_end:length(spk_x));
        spk_y_begin= spk_y(1:idx_begin);
        spk_y_end= spk_y(idx_end:length(spk_y));
        
        % get rate matrix - using function: Create_Rate_Map
        rate_mat_begin=CreateRateMap(pos_mean_x_begin,pos_mean_y_begin,pos_t_begin,spk_x_begin,spk_y_begin,spk_t_begin,parms);
        rate_mat_end=CreateRateMap(pos_mean_x_end,pos_mean_y_end,pos_t_end,spk_x_end,spk_y_end,spk_t_end,parms);
        rate_mat_total=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        %try
        
        % get zone mats of the half sessions
        % get max indices of full session
        [zone_mat_b, ~, ~]= getZones(rate_mat_begin);
        [zone_mat_e, ~, ~]= getZones(rate_mat_end);
        [~, peak_rates_t, max_inds]= getZones(rate_mat_total);
        
        %make sure full session max_inds not larger than half sessions
        size_b= size(zone_mat_b);
        size_e= size(zone_mat_e);
        if any(max_inds(:,1) > size_b(1))
            ind=find(max_inds(:,1) > size_b(1));
            max_inds(ind,1)= size_b(1);
        end
        
        if any(max_inds(:,1) > size_e(1))
            ind=find(max_inds(:,1) > size_e(1));
            max_inds(ind,1)= size_e(1);
        end
        
        if any(max_inds(:,2) > size_b(2))
            ind=find(max_inds(:,2) > size_b(2));
            max_inds(ind,2)= size_b(2);
        end
        
        if any(max_inds(:,2) > size_e(2))
            ind=find(max_inds(:,2) > size_e(2));
            max_inds(ind,2)= size_e(2);
        end
            
        % find peak rates of half sessions using full session max inds
        peak_rates_begin= nan(1,length(max_inds));
        for cen= 1:length(max_inds);
            peak_rates_begin(cen)= zone_mat_b(max_inds(cen,1), max_inds(cen,2));
        end
        
        peak_rates_end= nan(1,length(max_inds));
        for cen= 1:length(max_inds);
            peak_rates_end(cen)= zone_mat_e(max_inds(cen,1), max_inds(cen,2));
        end
        
        % put peak_rates in order of increasing order from total session
        ordered_inds= nan(1,length(peak_rates_t));
        sorted_peaks=sort(peak_rates_t);
        
        for h=1:length(peak_rates_t)
            ordered_inds(h)= find(sorted_peaks(h) == peak_rates_t);
        end
        
        peak_rates_b{i}=peak_rates_begin(ordered_inds);
        peak_rates_e{i}=peak_rates_end(ordered_inds);
        
        cd('N:\users\rebekkah\results and info of analysis')
        save('peak rates of half sessions', 'peak_rates_b', 'peak_rates_e');
        cd(parms.dir_load_data);
        
        
        %end %end of try
    end %end of if not fifth file (corrupted file)
end
