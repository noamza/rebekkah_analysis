parms.dir_load_data = 'C:\Users\Dori\Desktop\PF 8 results with wall scores'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('C:\Users\Dori\Desktop\Rebekkah data\module_list_distances_PF_8');

list=1;

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};
    dat = load(file_name);
    Cell1= dat.S;
    
    for l=k+1 : length(file_names)
        file_name2 = file_names{l};
        dat = load(file_name2);
        Cell2= dat.S;
        
        if Cell1.smallest_degree_wall == 1 | Cell1.smallest_degree_wall == 3
            smallest_degree_wall_1 = 1;
        elseif Cell1.smallest_degree_wall == 2 | Cell1.smallest_degree_wall == 4
            smallest_degree_wall_1 = 2;
        end
            
        if Cell2.smallest_degree_wall == 1 | Cell2.smallest_degree_wall == 3
            smallest_degree_wall_2 = 1;
        elseif Cell2.smallest_degree_wall == 2 | Cell1.smallest_degree_wall == 4
            smallest_degree_wall_2 = 2;
        end
                
        if smallest_degree_wall_1 == smallest_degree_wall_2
            orientation_diff= 0;
            angle_diff= abs(Cell1.smallest_degree - Cell2.smallest_degree);
        elseif smallest_degree_wall_1 ~= smallest_degree_wall_2
            orientation_diff= 10;
            angle_diff= NaN;
        end
           
            
   module_list(list,6) = orientation_diff;
   module_list(list,7) = angle_diff;
   
   list = list +1;
   
    
    end
end

save('module_orientation_list_PF8', 'module_list');


%% add 8th column same or not same based on distance (0.6)
for h= 1:length(module_list);
    if module_list(h,5) <= 0.5
        module_list (h, 8) = 2;
    elseif module_list(h,5) <0.7 & module_list(h,5) >0.5
        module_list(h,8) = 1;
    elseif module_list(h,5) > 0.7
        module_list (h,8) = 0;
    end
end

count1 =1;
count2=1;

for h= 1: length(module_list);
    if module_list(h,6)== 0
        same_orientation(count1) = module_list(h,4);
        same_orient_dist(count1) = module_list(h,8);
        
        count1=count1+1;
    elseif module_list(h,6)== 10
        diff_orientation(count2) = module_list(h,4);
        diff_orient_dist(count2) = module_list(h,8);
        
        count2=count2+1;
    end
end


figure;
n= 2
m= 2

subplot(n,m,1)
hist(same_orientation(:));

subplot(n,m,2)
hist(diff_orientation(:));

subplot(n,m,3)
hist(same_orient_dist(:));

subplot(n,m,4)
hist(diff_orient_dist(:));

disp('');
