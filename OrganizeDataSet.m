parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\Data - Original';
parms.dir_Sar='\\192.114.21.198\Dori_Data\data\rebekkah\original data- divided\Sargolini Data';
parms.dir_Der='\\192.114.21.198\Dori_Data\data\rebekkah\original data- divided\Derdikman Data';
parms.dir_Bon='\\192.114.21.198\Dori_Data\data\rebekkah\original data- divided\Bonnevie Data';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for i=453:length(file_names)
    file_name = file_names{i};
    dat=load(file_name);
    Cell=dat.S;
    
    if isfield(Cell, 'directory') 
        if findstr(Cell.directory,'Sargolini')
            copyfile(file_name, parms.dir_Sar);
            cd(parms.dir_load_data);
        elseif findstr(Cell.directory, 'bonnevie')
            copyfile(file_name, parms.dir_Bon);
            cd(parms.dir_load_data);
        end
    else
        copyfile(file_name, parms.dir_Der);
        cd(parms.dir_load_data);
    end
    
end