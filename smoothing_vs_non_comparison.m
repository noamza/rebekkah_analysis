parms.dir_load_data = 'C:\Users\Dori\Desktop\smoothing sigma 1';
parms.dir_save_pictures = 'C:\Users\Dori\Desktop\smoothing 1 photos';


    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    distances = [];
    
    for n=1:length(file_names)
        
        file_name = file_names{n};
        load(file_name);
        
       maxs(i,1) = find(mean_values_list==max(mean_values_list));  
       maxs(i,2) = find(peak_rates_list==max(peak_rates_list));
       maxs(i,3)= find(peak_values_list==max(peak_values_list));
       
       fig = figure;
       
        subplot(1,3,1);
        imagesc(old_zone_mat);
        str = sprintf('peak = %f',maxs(i,2));
        title(str);
        
        subplot(1,3,2);
       imagesc(new_zone_mat);
       str2 = sprintf('peak = %f',maxs(i,1));
       title(str2);
       
    subplot(1,3,3);
       imagesc(peak_zone_mat);
       str2 = sprintf('peak = %f',maxs(i,3));
       title(str2);
       
       
        cd(parms.dir_save_pictures);
      % saveas(fig,strcat('pic_',file_name, '.fig')); 
       saveas(fig,strcat('pic_',file_name, '.jpg')); %    
       cd(parms.dir_load_data);
       
       
    end
    
    disp('');