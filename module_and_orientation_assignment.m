parms.dir_load_data = 'C:\Users\Dori\Desktop\bin size 6 nonsmooth\results updated'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

modules=1;

for i =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{i};    
    dat = load(file_name); 
    Cell1= dat.S;

    for j=i+1 : length(file_names)
    file_name = file_names{j};    
    dat2 = load(file_name); 
    Cell2= dat2.S;
   
    if size(Cell1.rate_mat) == size(Cell2.rate_mat)
    
    
%% finds the modules difference (intersection/union)    
    
x10=Cell1.six_orientation_pts(7,1);
y10=Cell1.six_orientation_pts(7,2);
 beta = Cell1.module.phi ;
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha =0:pi/100:2*pi;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    x1 = x10 + (Cell1.module.major * cosalpha * cosbeta - Cell1.module.minor * sinalpha * sinbeta);
    y1 = y10 + (Cell1.module.major * cosalpha * sinbeta + Cell1.module.minor * sinalpha * cosbeta);

 x20=Cell2.six_orientation_pts(7,1);
y20=Cell2.six_orientation_pts(7,2);
 beta = Cell2.module.phi ;
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha =0:pi/100:2*pi;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    x2 = x20 + (Cell2.module.major * cosalpha * cosbeta - Cell2.module.minor * sinalpha * sinbeta);
    y2 = y20 + (Cell2.module.major * cosalpha * sinbeta + Cell2.module.minor * sinalpha * cosbeta);

if isreal(x1) && isreal(x2)&& isreal(y1) && isreal(y2)
[xa, ya] = polybool('union', x1, y1, x2, y2);
[xb, yb] = polybool('intersection', x1, y1, x2, y2);
 

 intersection = polyarea(xb,yb); 
 union=polyarea(xa,ya); 
else
   intersection = nan; 
 union=nan; 
end

module_diff = intersection/union; 

clear union;

module_list(modules, 1) = Cell1.i;
module_list(modules, 2) = Cell2.i;
module_list(modules, 3) = module_diff;


%%% finds if max peak is on the same or semi same wall location 

    location_inds1= [];
    location_inds2= [];
    
    location_inds1 = find(Cell1.location == 1);
    location_inds2= find(Cell2.location == 1);
    
    location1_len= length(location_inds1);
    location2_len= length(location_inds2);    
    
    if isempty(location_inds1) | isempty(location_inds2) 
        module_list(modules,4) = 5; %if one of the cells isn't touching border, don't count it
    elseif location1_len == location2_len & location_inds1 == location_inds2 
        module_list(modules,4) = 2; % if corner-corner or border-border same, same 
    elseif location1_len == 1 & location2_len== 1 & location_inds1 ~= location_inds2 
        module_list(modules,4) = 0; %if border-border and not same, different 
    elseif location1_len ==2 & location2_len ==2 & length(union(location_inds1,location_inds2)) == 3
        module_list(modules,4) = 0.5; %if corner-corner and touching one same border, semi-semi-same
    elseif location1_len ==2 & location2_len ==2 & length(union(location_inds1,location_inds2)) == 4
        module_list(modules,4) = 0; % if corner-corner and no same borders, different
    elseif length(union(location_inds1,location_inds2)) == 3
        module_list(modules,4) = 0; % if corner-border and no same borders, different
    elseif length(union(location_inds1,location_inds2)) == 2
        module_list(modules,4) = 1; %if corner-border and one same border, semi-same
    else
        disp('wtf');
    end    
    
    module_list(modules,5) = Distance(Cell1.norm_max_index(1,1), Cell1.norm_max_index(1,2), ...
            Cell2.norm_max_index(1,1), Cell2.norm_max_index(1,2)); 

        %% same location of max peak or not based on distance
        
     if module_list(modules,5) <= 0.5 
        module_list (modules,6) = 2;   %same
%     elseif module_list(modules,5) <0.7 & module_list(modules,5) >0.5
%         module_list(modules,6) = 0;
    elseif module_list(modules,5) > 0.5
        module_list (modules,6) = 0; %different
     end
     
    
     
   %% add 7 and 8 cloumn of module list. orientation wall same or diff and angle difference
     
         if Cell1.smallest_degree_Wall == 1 | Cell1.smallest_degree_Wall == 3
            smallest_degree_wall_1 = 1;
        elseif Cell1.smallest_degree_Wall == 2 | Cell1.smallest_degree_Wall == 4
            smallest_degree_wall_1 = 2;
        end
            
        if Cell2.smallest_degree_Wall == 1 | Cell2.smallest_degree_Wall == 3
            smallest_degree_wall_2 = 1;
        elseif Cell2.smallest_degree_Wall == 2 | Cell2.smallest_degree_Wall == 4
            smallest_degree_wall_2 = 2;
        end
                
        if smallest_degree_wall_1 == smallest_degree_wall_2
            orientation_diff= 0;
            angle_diff= abs(Cell1.smallest_degree - Cell2.smallest_degree);
        elseif smallest_degree_wall_1 ~= smallest_degree_wall_2
            orientation_diff= 10;
            angle_diff= NaN;
        end
           
            
   module_list(modules,7) = orientation_diff;
   module_list(modules,8) = angle_diff;
      
    modules = modules+1;

        end
    
    end
end 
save('module orientation list all cell pairs updated results', 'module_list');

figure; hist(module_list(:,3));     %all modules scores 
figure; hist(module_list(:,6));     %same or different max peak locations

 same_module_locations = [];
     non_same_module_locations=[];
     semi_same_module_locations = [];
     
    count=1;
    count2=1;
    count3=1;
    count4=1;
    
for modules = 1:length(module_list)     
     if module_list(modules,3) < 0.1
         non_same_module_locations(count) = module_list(modules,6);
         count = count+1;
         %     elseif module_list(k,3) > 0.6 && module_list(k,3) < 0.8
         %         semi_semi_same_module_locations(count4) = module_list(k,6);
         count4=count4+1;
     elseif module_list(modules,3) > 0.1 && module_list(modules,3) < 0.9
         semi_same_module_locations(count2) = module_list(modules,6);
         count2=count2+1;
     elseif module_list(modules,3) > 0.9
         same_module_locations(count3) = module_list(modules,6);
         count3=count3+1;       
     end
end

figure;

n= 1;
m= 3;

subplot(n,m,1)
%hist(same_module_locations);
bar(hist(same_module_locations) ./ sum(hist(same_module_locations)))
title('Same Module')
ylim([0 0.8])

subplot(n,m,3)
%hist(non_same_module_locations);
bar(hist(non_same_module_locations) ./ sum(hist(non_same_module_locations)))
title('Different Module')
ylim([0 0.8])
% 
subplot(n,m,2)
bar(hist(semi_same_module_locations) ./ sum(hist(semi_same_module_locations)))
title('Semi Same')
ylim([0 0.8])


remove = find(module_list(:,4)==5);
module_list(remove, 4) = NaN;

count1=1;
count2=1;
for modules = 1:length(module_list);
    if module_list(modules,7)== 0
        same_orientation(count1) = module_list(modules,4);
        same_orient_dist(count1) = module_list(modules,6);
        
        count1=count1+1;
    elseif module_list(modules,7)== 10
        diff_orientation(count2) = module_list(modules,4);
        diff_orient_dist(count2) = module_list(modules,6);
        
        count2=count2+1;
    end
end


figure;
n= 2;
m= 2;

subplot(n,m,1)
hist(same_orientation(:));

subplot(n,m,2)
hist(diff_orientation(:));

subplot(n,m,3)
hist(same_orient_dist(:));

subplot(n,m,4)
hist(diff_orient_dist(:));

count1=1;
count2=1;

for list = 1:length(module_list);
    if module_list(list, 8) <= 5
        same_angle(count1) = module_list(list,4);
        count1=count1+1;
    elseif module_list(list,8) >= 5
        diff_angle(count2) = module_list(list,4);
        count2=count2+1;
    end
end
        
count1=1;
count2=1;       
        
for list = 1:length(module_list);
    if module_list(list, 8) <= 5
      same_angle_dist(count1) = module_list(list,6);
          count1=count1+1;
    elseif module_list(list,8) >= 5
       diff_angle_dist(count2) = module_list(list,6);
       count2=count2+1;
    end
end

figure;

subplot(2,2,1)
hist(same_angle(:));

subplot(2,2,2)
hist(diff_angle(:));

subplot(2,2,3)
hist(same_angle_dist(:));

subplot(2,2,4)
hist(diff_angle_dist(:));

disp('')



count1=1;
count2=1;
count3=1;
count4=1;

for modules=1:length(module_list)
    if module_list(modules,3) > 0.7 & module_list(modules,8) < 2
    same_same(count1)=module_list(modules,6);
    same_same_dist(count1)=module_list(modules,5);
    count1=count1+1;
    elseif module_list(modules,3) < 0.7 & module_list(modules,8) < 2
    diff_same(count2)=module_list(modules,6);
    diff_same_dist(count2)=module_list(modules,5);
    count2=count2+1;
    elseif module_list(modules,3) > 0.7 & module_list(modules,8) > 2
    same_diff(count3)=module_list(modules,6);
    same_diff_dist(count3)=module_list(modules,5);
    count3=count3+1;
    elseif module_list(modules,3) < 0.7 & module_list(modules,8) > 2
    diff_diff(count4)=module_list(modules,6);
    diff_diff_dist(count4)=module_list(modules,5);
    count4=count4+1;
    end
end

values(1) = sum(same_same==2) ./ length(same_same)
values(2) = sum(same_diff==2) ./ length(same_diff)
values(3)= sum(diff_diff==2) ./ length (diff_diff)
figure; bar(values(:))