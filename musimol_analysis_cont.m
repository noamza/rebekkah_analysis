parms.dir_load_data = 'N:\users\rebekkah\musimol results criteria increase';

cd(parms.dir_load_data);

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\DB*.mat'));
file_names = {dir_list.name};

count=1;
count2=1;

for i=1:length(file_names)
    file_name = file_names{i};
    load(file_name);
    

    if mean_firing(4)/mean_firing(1) <= 0.5 && firing_rates(4)/firing_rates(1) <=0.8
      cluster_d(count)= cluster_score; 
      gridness_d(count)= gridness_scores(1)-gridness_scores(4);
      firing_rate_d(count)= firing_rates(4)/firing_rates(1);
       mean_firing_d(count)= mean_firing(4)/mean_firing(1);
      count=count+1;
    elseif mean_firing(4)/mean_firing(1) > 0.5 
     cluster_nd(count2)= cluster_score; 
      gridness_nd(count2)= gridness_scores(1)-gridness_scores(4);
      firing_rate_nd(count2)= firing_rates(4)/firing_rates(1);
       mean_firing_nd(count2)= mean_firing(4)/mean_firing(1);
       file_name
      count2=count2+1;
    else
        disp('huh')
        firing_rates(4)
    end
    
    total_firing(i)= firing_rates(4)/firing_rates(1);
    total_cluster(i)= cluster_score; 
    total_gridness(i)= gridness_scores(1)-gridness_scores(4);
    
end

mean(cluster_d)
mean(cluster_nd)

mean(gridness_d)
mean(gridness_nd)

mean(firing_rate_d)
mean(firing_rate_nd)

mean(mean_firing_d)
mean(mean_firing_nd)