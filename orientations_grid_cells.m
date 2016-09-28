parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results 202 cells';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
       
orientations(h) = S.smallest_degree;

orientations_improved(h)= S.min_angle;

end

figure; hist(orientations); 
figure; hist(orientations_improved); title('improved angle analysis')
