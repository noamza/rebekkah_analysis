% shuffle max peaks and then measure distance between pairs of cells.

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\FINAL ADJUSTED by rat and arena\17';

%load('N:\users\rebekkah\final data smoothed\FINAL ADJUSTED by rat and arena\module orientation same rat adjusted angles.mat')

load('N:\users\rebekkah\final data smoothed\FINAL ADJUSTED by rat and arena\16\shuffled module list.mat')

%list=1;

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

%%comment after first run

%  shuffled_module_list=module_list;
%  shuffled_module_list(:,9)= NaN;
%  shuffled_module_list(:,10)= NaN;

for i =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell1= dat.S;
    
    for j=i+1 : length(file_names)
        file_name = file_names{j};
        dat2 = load(file_name);
        Cell2= dat2.S;
        
        % shuffles max peak location- cell1
        
        [size_x, size_y] = size(Cell1.peak_zone_mat);
        
        peak_rates_list=[];
        
        peak_rates_list = Cell1.peak_rates(randperm(length(Cell1.peak_rates)));
        
        zone_num = find(peak_rates_list==max(peak_rates_list));
        zone_num = zone_num(1);
        
        [row, col] = find(Cell1.peak_zone_mat == zone_num);
        
        shuffle_index_cell1 = [];
        
        [size_x, size_y]=size(Cell1.rate_mat);
        
        shuffle_index_cell1(1) = row/size_x;
        shuffle_index_cell1(2) = col/size_y;
        
        % shuffles max peak location- cell2
        
        [size_x, size_y] = size(Cell2.peak_zone_mat);
        
        peak_rates_list=[];
        
        peak_rates_list = Cell2.peak_rates(randperm(length(Cell2.peak_rates)));
        
        zone_num = find(peak_rates_list==max(peak_rates_list));
        zone_num = zone_num(1);
        
        [row, col] = find(Cell2.peak_zone_mat == zone_num);
        
        shuffle_index_cell2 = [];
        
        [size_x, size_y]= size(Cell2.rate_mat);
        
        shuffle_index_cell2(1) = row/ size_x; %normalized index
        shuffle_index_cell2(2) = col/ size_y;
        
        %finds distance between max peaks
        
        shuffled_module_list(list,9) = Distance(shuffle_index_cell1(1), shuffle_index_cell1(2), ...
            shuffle_index_cell2(1), shuffle_index_cell2(2));
        
        
        % same location of max peak or not based on distance
        
        if shuffled_module_list(list,9) < 0.5
            shuffled_module_list (list,10) = 2; %same
        elseif shuffled_module_list(list,9) >= 0.5
            shuffled_module_list (list,10) = 0; %different
        end
        
        list=list+1;
        
    end
end

shuffled_module_list

save('shuffled module list', 'shuffled_module_list')  

