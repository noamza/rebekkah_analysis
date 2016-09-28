parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);
    Cell=dat.S;

    rate_mat_thresh= Cell.rate_mat;
    cut_off= max(max(rate_mat_thresh));
    rate_mat_thresh(rate_mat_thresh < 0.8 * cut_off)= 0;  
    
    figure; imagesc(rate_mat_thresh);
   
    disp('');
end