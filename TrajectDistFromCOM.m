
% finds the distance of the trajectory from the closest center of mass at
% each timestamp

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

parms.bin_size=3;


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    dist_from_cen= nan(1,length(pos_mean_x));
    
     %convert position coordinates to rate map coordinates
     [posi_x, posi_y]= ConvertCoordinates(Cell.rate_mat, parms.bin_size, pos_mean_x,pos_mean_y);
    
    for h= 1:length(pos_mean_x);
        
        if ~isnan(pos_mean_x(h)) & ~isnan(pos_mean_y(h))
            
           
            
            %find closest field COM
            closest_cen_coord= FindClosestCenCoords(Cell.max_inds, posi_x(h), posi_y(h));
            
            %calculate distance from trajectory to COM
            dist_from_cen(h)= Distance(closest_cen_coord(1), closest_cen_coord(2), posi_x(h), posi_y(h));
            
        else
            dist_from_cen(h)=nan;
        end
        
        
    end
    
   figure;
    plot(Cell.pos.t, dist_from_cen); hold on;
    ylabel('distance from field COM');
    xlabel('time');
    
    timestamp_inds= nan(1,length(Cell.spk.t));
    
    for h=1:length(Cell.spk.t)
    inds= find(Cell.pos.t >= Cell.spk.t(h));
    timestamp_inds(h) = inds(1);
    end
    
    dist_at_spk= dist_from_cen(timestamp_inds);
    
    %plot time and distance from COM of all spikes
    plot(Cell.spk.t, dist_at_spk,  'r.', 'MarkerSize', 5); 
    
    
    disp('');
end