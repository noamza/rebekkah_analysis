parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

firing_vs_time_corr= nan(1,length(file_names));

parms.bin_size=6;

cd('N:\users\rebekkah\results and info of analysis')
load('variability and border distances UPDATED.mat')
cd(parms.dir_load_data)

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    % find the mean position of the two LEDs
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    % info for zone mat
    max_inds= max_indices{i};
    peak_rates= peak_rates_all{i};
    num_zone_mat= number_zone_mat{i};
        
    [pos_x_inds, pos_y_inds]= ConvertCoordinates(num_zone_mat, parms.bin_size, pos_mean_x,pos_mean_y) ;
    [spk_x_inds, spk_y_inds]= ConvertCoordinates(num_zone_mat, parms.bin_size, spk_x,spk_y);
    
    % find number of spikes per zone
    num_of_fields= length(max_inds);
    number_spks_per_field=zeros(1,num_of_fields);
    for h=1:length(spk_x_inds);
        if ~isnan(spk_x_inds(h))    % skip cases with nan
            zone_num= num_zone_mat(spk_x_inds(h), spk_y_inds(h));
            if zone_num ~= 0
                number_spks_per_field(zone_num)= number_spks_per_field(zone_num)+1;
            end
        end
    end
    
    %find time spent per zone
    time_spent_per_field=zeros(1,num_of_fields);
    for h=1:length(pos_x_inds);
        if ~isnan(pos_x_inds(h))    % skip cases with nan
            zone_num= num_zone_mat(pos_x_inds(h), pos_y_inds(h));
            if zone_num ~= 0
                time_spent_per_field(zone_num)= time_spent_per_field(zone_num)+1;
            end
        end
    end
    
    corr= corrcoef(peak_rates, time_spent_per_field);
    firing_vs_time_corr(i)= corr(2);
    
end

mean(firing_vs_time_corr)
median(firing_vs_time_corr)

figure; hist(firing_vs_time_corr);

disp('')