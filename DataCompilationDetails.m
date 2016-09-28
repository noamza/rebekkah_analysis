dbstop if error

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\data set for adaptation analysis- 0.4 grid 0.2 hd with circular';


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
   
  
    Bonnevie= 0;
    Derdikman=0;
    Sargolini=0;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);
    Cell=dat.S;
     
    
    if isfield(Cell, 'file_name') & strfind(Cell.file_name, 'DB_MUSC_MEC')
       
            Bonnevie= Bonnevie +1;
        
    elseif isfield(Cell, 'directory')
            if strfind(Cell.directory,  'Sargolini 2006')
                Sargolini= Sargolini +1;
            end
    else
        Derdikman= Derdikman +1;
    end
    
end

Bonnevie
Derdikman
Sargolini  