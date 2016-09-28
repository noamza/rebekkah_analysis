parms.dir_load_data = 'C:\Users\Dori\Desktop\gridness 0.3 HD 0.3 paraments'; 

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};    
    load(file_name); 
    
    orientations(i) = S.smallest_degree;
    
end

disp('');