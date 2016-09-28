
parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\rotated arenas by rat and arena AND DATE\19';

load('N:\users\rebekkah\final data smoothed\rotated arenas by rat and arena AND DATE\module orientation SAME DATE.mat')

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};    
    dat = load(file_name); 
    Cell1= dat.S;

    for j=k+1 : length(file_names)
    file_name = file_names{j};    
    dat2 = load(file_name); 
    Cell2= dat2.S;
    
    for modules= 1:length(module_list)
        if Cell1.i == module_list(modules,1) & Cell2.i == module_list(modules,2)
            %if module_list(modules,8) < 1 & module_list(modules,8) > 0
                 if Cell1.smallest_degree <= Cell2.smallest_degree
                     plot(Cell1.smallest_degree, Cell2.smallest_degree, 'x', 'markersize', 15)
                     hold on;
                 elseif Cell2.smallest_degree < Cell1.smallest_degree
                     plot(Cell2.smallest_degree, Cell1.smallest_degree, 'x', 'markersize', 15)
                     hold on;
                 end
            %end
        end
    end
    
    end
end