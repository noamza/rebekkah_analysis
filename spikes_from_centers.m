
parms.dir_load_data= 'N:\users\rebekkah\are spks closer to max peaks\results'; 



dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count=1;
count2=1;
% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    
    if peak_score >=1.5
    max_zone_distance(count)= ordered_distances(end);
    other_zones_mean_distance(count)= mean(ordered_distances(1:end-1));
    
    var_max_zone_distance(count)= std_ordered_distances(end);
    var_other_zones_mean_distance(count)= mean(std_ordered_distances(1:end-1));
    count=count+1;
    end


    if peak_score >=1.75
    max_zone_distance2(count2)= ordered_distances(end);
    other_zones_mean_distance2(count2)= mean(ordered_distances(1:end-1));
    
    var_max_zone_distance2(count2)= std_ordered_distances(end);
    var_other_zones_mean_distance2(count2)= mean(std_ordered_distances(1:end-1));
    count2=count2+1;
    end
    
end

x=0:6;
y=0:6;

figure; scatter(max_zone_distance, other_zones_mean_distance); hold on;
plot(x,y,'-')
xlabel('mean spike distance from maximum-firing field center', 'fontname', 'calibir', 'fontsize', 14)
ylabel('mean spike distance from non-max field centers', 'fontname', 'calibir', 'fontsize', 14)
figure; scatter(max_zone_distance2, other_zones_mean_distance2); hold on;
xlabel('mean spike distance from maximum-firing field center', 'fontname', 'calibir', 'fontsize', 14)
ylabel('mean spike distance from non-max field centers', 'fontname', 'calibir', 'fontsize', 14)
plot(x,y,'-')

