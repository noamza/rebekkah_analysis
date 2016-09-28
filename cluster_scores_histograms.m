parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS/third';


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);

cluster_scores_three(i)= cluster_score;

end

% fig=figure;
% hist(cluster_scores_two); 
% h = findobj(gca,'Type','patch');
% set(h, 'facecolor', 'r', 'facealpha', 0.75);

hold on;

.........second histogram

parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\second';


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);

cluster_scores_two(i)= cluster_score;

end



% hist(cluster_scores);
% h = findobj(gca,'Type','patch');
% set(h, 'facealpha', 0.75);

%save('cluster scores at least two arena', 'cluster_scores')


parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\max';


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);

cluster_scores(i)= cluster_score;

end


clusters= cluster_scores;
clusters(end+1:76)= cluster_scores_two;
clusters(end+1:end+38)= cluster_scores_three;

ordered(1:38)=1;
ordered(39:76)=2;
ordered(77:114)=3;

figure; boxplot(clusters, ordered);