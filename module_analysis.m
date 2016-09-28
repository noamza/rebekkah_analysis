

parms.dir_load_data = 'C:\Users\Dori\Desktop\elipse plots PF 8 results'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('C:\Users\Dori\Desktop\Rebekkah data\module list PF 8 half smooth');

modules =1;

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};    
    dat = load(file_name); 
    Cell1= dat.S;

    for l=k+1 : length(file_names)
    file_name2 = file_names{l};    
    dat = load(file_name2); 
    Cell2= dat.S;
   
    location_inds1= [];
    location_inds2= [];
    
    location_inds1 = find(Cell1.location == 1);
    location_inds2= find(Cell2.location == 1);
    
    location1_len= length(location_inds1);
    location2_len= length(location_inds2);
    
    if location1_len == 0 | location2_len == 0 
        module_list(modules,4) = 5; %if one of the cells isn't touching border, don't count it
    elseif location1_len == location2_len & location_inds1 == location_inds2 
        module_list(modules, 4) = 2; % if corner-corner or border-border same, same 
    elseif location1_len == 1 & location2_len== 1 & location_inds1 ~= location_inds2 
        module_list(modules,4) = 0; %if border-border and not same, different 
    elseif location1_len ==2 & location2_len ==2 & length(union(location_inds1,location_inds2)) == 3
        module_list(modules,4) = 0.5; %if corner-corner and touching one same border, semi-semi-same
    elseif location1_len ==2 & location2_len ==2 & length(union(location_inds1,location_inds2)) == 4
        module_list(modules,4) = 0; % if corner-corner and no same borders, different
    elseif length(union(location_inds1, location_inds2)) == 3
        module_list(modules,4) = 0; % if corner-border and no same borders, different
    elseif length(union(location_inds1, location_inds2)) == 2
        module_list(modules,4) = 1; %if corner-border and one same border, semi-same
    else
        disp('wtf');
    end    
    
    module_list(modules,5) = Distance(Cell1.norm_max_index(1,1), Cell1.norm_max_index(1,2), ...
            Cell2.norm_max_index(1,1), Cell2.norm_max_index(1,2)); 
    
    
    modules= modules+1;
    
    end
end

save('module_list_distances_PF_8_halfsmooth', 'module_list');

disp('');