
parms.dir_load_data = 'C:\Users\Dori\Desktop\results with diff autocorr\results';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
    orientations(h)= S.smallest_degree;
    
end

hist(orientations(:));

disp('');