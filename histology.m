

dbstop if error

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\combined 3 sets MD';
%parms.dir_save_pictures= 'N:\users\rebekkah\results and info of analysis\final images';

cd(parms.dir_load_data);

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for i =1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);

    if isfield(dat,'db') %Bonnevie
        histology{i}=dat.db.area; 
    elseif isfield('histology_layer',Cell) 
        histology{i}=Cell.histology_layer; 
    end
    
end