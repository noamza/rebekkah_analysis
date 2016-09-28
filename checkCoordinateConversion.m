parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

cd(parms.dir_load_data);

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};



for i=1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    max_inds_c=[];
    pos_mean_x=[];
    pos_mean_y=[];
    
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    
%convert max_inds to spike_mat coordinates

max_inds_c(:,1)= Cell.max_inds(:,1)*3 - max(pos_mean_x);
max_inds_c(:,2)= Cell.max_inds(:,2)*3 - max(pos_mean_y);
    
    figure;
    
    subplot(1,3,1)
    plot(spk_x, spk_y, '.'); hold on;
    plot(max_inds_c(:,2), max_inds_c(:,1), 'ok', 'MarkerFaceColor', 'k')
   
    
    subplot(1,3,2)
     plot(spk_x, spk_y, '.'); hold on;
    plot(max_inds_c(:,1), max_inds_c(:,2), 'ok', 'MarkerFaceColor', 'k')
    
    subplot(1,3,3)
    imagesc(Cell.rate_mat); hold on;
    plot(Cell.max_inds(:,2), Cell.max_inds(:,1), 'ok', 'MarkerFaceColor', 'k')
    
    disp('');
    
end