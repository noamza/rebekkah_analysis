parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
parms.dir_save_data = 'C:\Users\Dori\Desktop\only peak rates for shuffling';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

parms.bin_size= 3;
parms.sigma= 1.5;

for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    %    clear dat;
    %save speed_thresh_hold vector
    % Cell.speed_thresh_hold = speed_thresh_hold;
    
    if i ~= 5
        
        
        pos_t_len= length(Cell.pos.t);
        
        pos_t_begin= Cell.pos.t(1:15000);   %first half of session
        pos_t_end= Cell.pos.t(pos_t_len-14999:pos_t_len); %second half of session
        
        tmp_begin = abs(Cell.spk.t-pos_t_begin(15000)); %subtract last timestamp of first session
        tmp_end= abs(Cell.spk.t- pos_t_end(1)) ;% subtract first timestamp of second session
        [~, idx_begin] = min(tmp_begin); %index of first session
        [~, idx_end] = min(tmp_end) ;%index of first session
        spk_t_begin = Cell.spk.t(1:idx_begin); %spike vector of first sesh
        spk_t_end = Cell.spk.t(idx_end:length(Cell.spk.t)); %spike vector of second sesh
        
        if ~isempty(Cell.pos.x2)
            dt = Cell.pos.t(2)-Cell.pos.t(1);
            
            % calculate the the rat's head direction (using average of both leds)
            pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
            pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
            
            % build the axis of location when spikes are made
            spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
            spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
            
            % divide into first 5 minutes and last 5 minutes
            pos_mean_x_begin= pos_mean_x(1:15000); %x position for first half of session
            pos_mean_x_end= pos_mean_x(pos_t_len-14999:pos_t_len); %x position for second half of session
            pos_mean_y_begin= pos_mean_y(1:15000); %y position for first half of session
            pos_mean_y_end= pos_mean_y(pos_t_len-14999:pos_t_len); %y position for second half of session
            
            % confirm the Cell.spk.t and spk_x are same size and correlate to
            % eachother!!!!!!!!!!!!!!!
            
            spk_x_begin= spk_x(1:idx_begin);
            spk_x_end= spk_x(idx_end:length(spk_x));
            spk_y_begin= spk_y(1:idx_begin);
            spk_y_end= spk_y(idx_end:length(spk_y));
            
            % get rate matrix - using function: Creat_Rate_Map
            rate_mat_begin=CreateRateMap(pos_mean_x_begin,pos_mean_y_begin,pos_t_begin,spk_x_begin,spk_y_begin,spk_t_begin,parms);
            rate_mat_end=CreateRateMap(pos_mean_x_end,pos_mean_y_end,pos_t_end,spk_x_end,spk_y_end,spk_t_end,parms);
            rate_mat_total=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
            
            S_begin=get_zones(rate_mat_begin, dat.S);
            S_end=get_zones(rate_mat_end, dat.S);
            S_total=get_zones(rate_mat_total, dat.S);
            
            peak_rates_begin=nan(1,length(Cell.peak_rates));
           peak_rates_end=nan(1,length(Cell.peak_rates));
            
            for h= 1:length(Cell.peak_rates)
                [ind1,ind2]=find(Cell.peak_zone_mat==h);
                
                [sizex,sizey]= size(S_begin.zone_mat);
                if sizey < ind2;
                    ind2=sizey;
                end
                if sizex < ind1;
                    ind1=sizex;
                end
                
                [sizex,sizey]= size(S_end.zone_mat);
                if sizey < ind2;
                    ind2=sizey;
                end
                if sizex < ind1;
                    ind1=sizex;
                end
                
                peak_rates_begin(h)=S_begin.zone_mat(ind1,ind2);
                peak_rates_end(h)= S_end.zone_mat(ind1,ind2);
            end
            
            ordered=1:length(Cell.peak_rates);
            ordered_inds= nan(1,length(Cell.peak_rates));
            
            for h=1:length(Cell.sorted_means)
                ordered_inds(h)= find(Cell.sorted_means(h) == Cell.peak_rates);
            end
            
            peak_rates_begin=peak_rates_begin(ordered_inds);
            peak_rates_end=peak_rates_end(ordered_inds);
            
            S= Cell;
            
            cd(parms.dir_save_data);
            save(sprintf('%s', file_name), 'peak_rates_begin', 'peak_rates_end');
            cd(parms.dir_load_data);
            
        end
    end
end
