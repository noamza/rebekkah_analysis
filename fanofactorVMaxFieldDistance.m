parms.dir_load_data = 'N:\users\rebekkah\bin size 6 nonsmooth\results updated';


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
    
    fano_factor(i)= var(Cell.sorted_means)/mean(Cell.sorted_means);
    
    [size_x,~]= size(Cell.rate_mat);
    
    distance(i)= Cell.max_peak_distance/size_x;
    
end

figure; scatter(distance, fano_factor, 'd');
xlabel('normalized distance', 'fontsize', 14, 'fontname', 'calibri')
xlabel('Fano factor', 'fontsize', 14, 'fontname', 'calibri')
hold on;
lsline