parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\FINAL ADJUSTED by rat and arena\17';  

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('N:\users\rebekkah\final data smoothed\FINAL ADJUSTED by rat and arena\module orientation same rat adjusted angles.mat');

count =count;

%count=1;

%figure;

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};    
    dat = load(file_name); 
    Cell1= dat.S;

    for l= k+1 : length(file_names)
    file_name2 = file_names{l};    
    dat = load(file_name2); 
    Cell2= dat.S;
    
    if module_list(count,3) >0.7
        if Cell1.smallest_degree_adj > Cell2.smallest_degree_adj
            scatter(Cell1.smallest_degree_adj, Cell2.smallest_degree_adj, 'markerfacecolor', 'r', 'markeredgecolor', 'r')
            hold on;
        elseif Cell2.smallest_degree_adj > Cell1.smallest_degree_adj
            scatter(Cell2.smallest_degree_adj, Cell1.smallest_degree_adj, 'markerfacecolor', 'r', 'markeredgecolor', 'r')
        hold on;
        end
    elseif module_list(count,3) <0.7
        if Cell1.smallest_degree_adj > Cell2.smallest_degree_adj
            scatter(Cell1.smallest_degree_adj, Cell2.smallest_degree_adj, 'markerfacecolor', 'b', 'markeredgecolor', 'b')
            hold on;
        elseif Cell2.smallest_degree_adj > Cell1.smallest_degree_adj
        scatter(Cell2.smallest_degree_adj, Cell1.smallest_degree_adj, 'markerfacecolor', 'b', 'markeredgecolor', 'b')
        end
    end
        
    count=count+1;
    
    end
end