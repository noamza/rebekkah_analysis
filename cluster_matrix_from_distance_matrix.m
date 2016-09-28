%create cluster matrix using distance matrix

parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\third max results above 0.3';
parms.dir_save_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\third cluster scores'; 

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);

    anchor_pt_distances=third_anchor_pt_distances;
    

   distance_matrix = anchor_pt_distances;
   

    cluster_matrix= nan(5);

    for g= 1:length(distance_matrix)-1;
        for h= g+1:length(distance_matrix);
            if distance_matrix(g,h) <= 0.35
                cluster_matrix(g,h)= 1;
            elseif distance_matrix(g,h) > 0.35
                cluster_matrix(g,h)= 0;
            end
        end
    end

cluster_score= sum(sum(cluster_matrix==1))/ sum(sum(~isnan(cluster_matrix)));

     cd(parms.dir_save_data);
            save(sprintf('%s.mat', file_name), 'cluster_matrix', 'cluster_score');
            cd(parms.dir_load_data);



end
