parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
parms.dir_save_pictures= 'N:\users\rebekkah\results and info of analysis\final images';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

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
    
    rate_mat=all_rate_mats{i};
    number_zone_mat=number_zone_mat_sm{i};
    
    [pos_x_inds, pos_y_inds]= ConvertCoordinates(rate_mat, 3, pos_mean_x,pos_mean_y);
    [spk_x_inds, spk_y_inds]= ConvertCoordinates(rate_mat, 3, spk_x,spk_y);
    field_t = FindFieldPerTimeBin(Cell.pos_t, pos_x_inds, pos_y_inds, number_zone_mat);
    
    num_of_fields= length(peak_rates_all_sm{i});
    number_spks_per_field=zeros(1,num_of_fields);
    for h=1:length(spk_x);
        zone_num= number_zone_mat 
        
    end
    
    
end