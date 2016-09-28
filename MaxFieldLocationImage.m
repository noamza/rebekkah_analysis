parms.dir_load_data = 'N:\users\rebekkah\bin size 6 nonsmooth\results updated';

%parms.bin_size=3;

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
    
    total_max_inds(i,1)= Cell.norm_max_index(1);  
    total_max_inds(i,2)= Cell.norm_max_index(2);
    
end

[max_pt_rate_map]= CreateMaxFieldRateMap(total_max_inds);

figure; imagesc(max_pt_rate_map);

sum(sum(max_pt_rate_map))

parms.sigma= 1.5;

max_pt_rate_map_smooth= SmoothRateMat(max_pt_rate_map, parms);

figure; imagesc(max_pt_rate_map_smooth);