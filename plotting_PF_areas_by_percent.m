
parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results 202 cells with adjusted zone mats [use to find PF sizes] [10 percent removed]';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

figure;

all_PF_areas_by_percent=[];

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};
    load(file_name)
    
    PF_areas_by_percent= S.PF_areas/max(S.PF_areas);
    
    if PF_areas_by_percent(end) < 0.7
    plot(S.norm_means, PF_areas_by_percent, 'o', 'color', 'r')
    hold on;
    else
      plot(S.norm_means, PF_areas_by_percent, 'o', 'color', 'b')
    hold on;
    end
    
    all_PF_areas_by_percent= [all_PF_areas_by_percent, PF_areas_by_percent];  
    
    max_peaks(k)= S.PF_areas(end)/max(S.PF_areas);
    
end

figure; hist(all_PF_areas_by_percent);