parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\202 cells by rat and arena\22'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('N:\users\rebekkah\final data smoothed\data sets\202 cells by rat and arena\21\module orientation SAME DATE.mat')
% 
 [sizex, sizey] = size(module_list);
% % 
  modules = sizex +1;
%
%modules=1;

for i =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{i};    
    dat = load(file_name); 
    Cell1= dat.S;

    for j=i+1 : length(file_names)
    file_name = file_names{j};    
    dat2 = load(file_name); 
    Cell2= dat2.S;
        
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
        module_list(modules, 4) = 2; % if corner-corner or border-border same, same 
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
    
    %%finds distance between max peaks
    
    module_list(modules,5) = Distance(Cell1.norm_max_index(1), Cell1.norm_max_index(2), ...
            Cell2.norm_max_index(1), Cell2.norm_max_index(2)); 

 
        
        %% same location of max peak or not based on distance
        
     if module_list(modules,5) < 0.5 
        module_list (modules,6) = 2; %same
%      elseif module_list(modules,5) > 0.4 & module_list(modules,5) < 0.6 %semisame
%          module_list(modules,6) = 0;
    elseif module_list(modules,5) >= 0.5
        module_list (modules,6) = 0; %different
     end
     
    
     
   %% add 7 and 8 cloumn of module list. orientation wall same or diff and angle difference
     
   
   
        if Cell1.smallest_degree_wall == 1 | Cell1.smallest_degree_wall == 2
            smallest_degree_wall_1 = 1;
        elseif Cell1.smallest_degree_wall == 3 | Cell1.smallest_degree_wall == 4
            smallest_degree_wall_1 = 2;
        end
            
        if Cell2.smallest_degree_wall == 1 | Cell2.smallest_degree_wall == 2
            smallest_degree_wall_2 = 1;
        elseif Cell2.smallest_degree_wall == 3 | Cell2.smallest_degree_wall == 4
            smallest_degree_wall_2 = 2;
        end
       
        % finds angle difference if angle is same orientation (top-bottom,
        % or right-left)
        
        if smallest_degree_wall_1 == smallest_degree_wall_2
            orientation_diff= 0;
            angle_diff= abs(Cell1.smallest_degree - Cell2.smallest_degree);
        elseif smallest_degree_wall_1 ~= smallest_degree_wall_2
            orientation_diff= 10;
            angle_diff= NaN;
        end
           
            
   module_list(modules,7) = orientation_diff;
   module_list(modules,8) = angle_diff;
   
 %  module_list(modules,9) = abs(Cell1.major_angle-Cell2.major_angle);
   
    modules = modules+1;

    end
    
end

save('module orientation SAME DATE', 'module_list');

module_list 

clearvars

%      same_module_locations = [];
%      non_same_module_locations=[];
%      semi_same_module_locations = [];
%      
%      
%     count=1;
%     count2=1;
%     count3=1;
%     count4=1;
%     count5=1;
%     
% for modules = 1:length(module_list)     
%      if module_list(modules,3) < 0.6
%          non_same_module_locations(count) = module_list(modules,6);
%          count = count+1;
%   %          elseif module_list(modules,3) > 0.7 && module_list(modules,3) < 0.8
%    %              semi_semi_same_module_locations(count4) = module_list(modules,6);
%  %        count4=count4+1;
%   %   elseif module_list(modules,3) > 0.8 && module_list(modules,3) < 0.9
%    %       middle_module_locations(count5)= module_list(modules,6);
%    %       count5=count5+1; 
%       elseif module_list(modules,3) > 0.6 && module_list(modules,3) < 0.9
%           semi_same_module_locations(count2) = module_list(modules,6);
%            count2=count2+1;
%      elseif module_list(modules,3) > 0.9
%          same_module_locations(count3) = module_list(modules,6);
%          count3=count3+1;       
%      end
% end
% 
% module_location =[];
% 
% module_location(1) = sum(non_same_module_locations==2)/length(non_same_module_locations);
% %module_location(2) = sum(semi_semi_same_module_locations==2)/length(semi_semi_same_module_locations);
% %module_location(3) = sum(middle_module_locations==2)/length(middle_module_locations);
% module_location(2) = sum(semi_same_module_locations==2)/length(semi_same_module_locations);
% module_location(3) = sum(same_module_locations==2)/length(same_module_locations);
% 
% 
% figure; bar(module_location);
% 
% figure;
% 
% n= 1;
% m= 3;
% 
% subplot(n,m,1)
% %hist(same_module_locations);
% bar(hist(same_module_locations) ./ sum(hist(same_module_locations)))
% title('Same Module')
% ylim([0 0.8])
% 
% subplot(n,m,3)
% %hist(non_same_module_locations);
% bar(hist(non_same_module_locations) ./ sum(hist(non_same_module_locations)))
% title('Different Module')
% ylim([0 0.8])
% % 
% subplot(n,m,2)
% bar(hist(semi_same_module_locations) ./ sum(hist(semi_same_module_locations)))
% title('Semi Same')
% ylim([0 0.8])
% % 
% % subplot(n,m,3)
% % bar(hist(semi_semi_same_module_locations) ./ sum(hist(semi_semi_same_module_locations)))
% % title('Semi Smi Same')
% % ylim([0 0.8])
% 
% remove = find(module_list(:,4)==5);
% module_list(remove, 4) = NaN;
% 
% change=find(module_list(:,4)==0.5 | module_list(:,4)== 1);
% module_list(change,4) = 2;
% 
% count1=1;
% count2=1;
% for modules = 1:length(module_list);
%     if module_list(modules,9)<= 5
%         %same_orientation(count1) = module_list(modules,4);
%         same_orient_dist(count1) = module_list(modules,6);
%         
%         count1=count1+1;
%     elseif module_list(modules,9)>5
%         %diff_orientation(count2) = module_list(modules,4);
%         diff_orient_dist(count2) = module_list(modules,6);
%         
%         count2=count2+1;
%     end
% end
% 
% figure;
% n= 2;
% m= 1;
% 
% % subplot(n,m,1)
% % hist(same_orientation(:));
% % 
% % subplot(n,m,2)
% % hist(diff_orientation(:));
% % 
% % subplot(n,m,3)
% % hist(same_orient_dist(:));
% % 
% % subplot(n,m,4)
% % hist(diff_orient_dist(:));
% 
% orient_loc_dist(1) = sum(same_orient_dist==2)/length(same_orient_dist);
% orient_loc_dist(2)= sum(diff_orient_dist==2)/length(diff_orient_dist);
% 
% figure; bar(orient_loc_dist)
% 
% 
% count1=1;
% count2=1;
% 
% for list = 1:length(module_list);
%     if module_list(list, 8) <= 5
%         same_angle(count1) = module_list(list,4);
%         count1=count1+1;
%     elseif module_list(list,8) > 5
%         diff_angle(count2) = module_list(list,4);
%         count2=count2+1;
%     end
% end
%         
% count1=1;
% count2=1;       
%         
% for list = 1:length(module_list);
%     if module_list(list, 8) <= 5
%       same_angle_dist(count1) = module_list(list,6);
%           count1=count1+1;
%     elseif module_list(list,8) > 5
%        diff_angle_dist(count2) = module_list(list,6);
%        count2=count2+1;
%     end
% end
% 
% figure;
% 
% subplot(2,2,1)
% hist(same_angle(:));
% 
% subplot(2,2,2)
% hist(diff_angle(:));
% 
% subplot(2,2,3)
% hist(same_angle_dist(:));
% 
% subplot(2,2,4)
% hist(diff_angle_dist(:));
% 
% 
% %% both same module and same orientation, same m.p. location?
% count1=1;
% count2=1;
% for modules = 1:length(module_list)
% if module_list(modules, 3)>= 0.7 & module_list(modules,8) <= 5
%     same_parameters_dist(count1)= module_list(modules, 5);
%    same_parameters(count1) = module_list(modules,6);
%     count1=count1+1;
% else
%      diff_parameters_dist(count2)= module_list(modules, 5);
%    diff_parameters(count2) = module_list(modules,6);
%     count2=count2+1;
% end
% end
% 
% count1=1;
% count2=1;
% for modules= 1:length(module_list)
%     if module_list(modules,3) >=0.8
%         same_modules_dist(count1)=module_list(modules,5);
%         count1=count1+1;
%     elseif module_list(modules,3) < 0.8 
%         diff_modules_dist(count2)=module_list(modules,5);
%         count2=count2+1;
%     end
% end
% 
% figure; hist(same_modules_dist); hold on;
% hist(diff_modules_dist);
% h = findobj(gca,'Type','patch');
% set(h(2),'Facecolor','r','EdgeColor','k');
% 
% count1=1;
% count2=1;
% count3=1;
% 
% for modules=1:length(module_list)
%     if module_list(modules,8) <2
%         same_orientation_distances(count1)= module_list(modules,5);
%         count1=count1+1;
%     elseif module_list(modules,8) >2 & module_list(modules,8)< 5
%         similar_orientation_distances(count2)=module_list(modules,5);
%         count2=count2+1;
%     elseif module_list(modules,8) >5
%         diff_orientation_distances(count3)=module_list(modules,5);
%         count3=count3+1;
%     elseif isnan(module_list(modules,8))
%         disp('NaN');
%     else
%         disp('wtf');
%     end
% end
% 
% figure; subplot(1,3,1);
% bar(hist(same_orientation_distances) ./ sum(hist(same_orientation_distances)))
% % subplot(1,3,2)
% % hist(similar_orientation_distances);
% subplot(1,3,3)
% bar(hist(diff_orientation_distances) ./ sum(hist(diff_orientation_distances)))
% 
% 
% 
% 
% for modules=1:length(module_list)
%     if module_list(modules,5) <= 0.3 
%         module_list (modules,6) = 2; %same
%      elseif module_list(modules,5) > 0.3 & module_list(modules,5) < 0.7 %semisame
%          module_list(modules,6) = 2;
%     elseif module_list(modules,5) > 0.7 & module_list(modules,5) <0.9
%         module_list(modules,6) = 0;
%     elseif module_list(modules,5) >= 0.9
%         module_list (modules,6) = 0; %different
%      end
% end
% 
% %same or not max peak whether same or not angle
% 
% count1=1;
% count2=1;
% count3=1;
% for modules=1:length(module_list)
%     if module_list(modules,8) <2
%         same_orientation_distances(count1)= module_list(modules,6);
%         count1=count1+1;
% %     elseif module_list(modules,8) >2 & module_list(modules,8)< 5
% %         similar_orientation_distances(count2)=module_list(modules,6);
% %         count2=count2+1;
%     elseif module_list(modules,8) >2
%         diff_orientation_distances(count3)=module_list(modules,6);
%         count3=count3+1;
%     elseif isnan(module_list(modules,8))
%         disp('NaN');
%     else
%         disp('wtf');
%     end
% end
% 
% display(1)= sum(same_orientation_distances==2)./sum(hist(same_orientation_distances))
% display(2)= sum(diff_orientation_distances==2)./sum(hist(diff_orientation_distances))
% figure; bar(display);
% 
% count1=1;
% count2=1;
% count3=1;
% count4=1;
% for modules=1:length(module_list)
%     if module_list(modules,3) >= 0.7 & module_list(modules,8) <= 2 %same module and same angle
%         all_same(count1) = module_list(modules,6);
%         all_same_dist(count1) = module_list(modules,5);
%         count1=count1+1;
%     elseif module_list(modules,3) >=0.7 & module_list(modules,8)>2 %same module, diff angle
%         module_same_angle_diff(count2) =module_list(modules,6);
%         module_same_angle_diff_dist(count2) =module_list(modules,5);
%         count2=count2+1;
%     elseif module_list(modules,3) <0.7 & module_list(modules,8)>2 %diff module and angle
%         all_diff(count3) =module_list(modules,6);
%         all_diff_dist(count3) =module_list(modules,5);
%         count3=count3+1;;
%     elseif module_list(modules,3) <0.7 & module_list(modules,8)<=2 %diff module, same angle
%         module_diff_angle_same(count4) =module_list(modules,6);
%         module_diff_angle_same_dist(count4) =module_list(modules,5);
%         count4=count4+1;
%     else
%         disp('huh')
%     end
% end
%  
% 
% %% comparing diff module and same orient, diff module and diff orient, and same module and same orient
% 
% count1=1;
% count2=1;
% count3=1;
% 
% for modules=1:length(module_list)
%     if module_list(modules,3) > 0.7 & module_list (modules,8) >4
%         same_m_diff_o(count1)= module_list(modules,6);
%         count1=count1+1;
%     elseif module_list(modules,3) <0.7 & module_list(modules,8) >4
%         diff_m_diff_o(count2)=module_list(modules,6)
%         count2=count2+1;
%     elseif module_list(modules,3) >0.7 & module_list(modules,8) <=4
%         same_m_same_o(count3)=module_list(modules,6)
%         count3=count3+1;
%     end
% end
% 
% orient_module(1) = sum(diff_m_diff_o ==2)/length(diff_m_diff_o)
% orient_module(2) = sum(same_m_diff_o ==2)/length(same_m_diff_o)
% orient_module(3) = sum(same_m_same_o ==2)/length(same_m_same_o)
% 
% figure; bar(orient_module);
% 
% %plot scatterplot of module diff vs angle diff by color depending on same
% %or diff modules
% 
% 
% figure;
% 
% for modules=1:length(module_list);
%     if module_list(modules,5) <0.4
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'r')
%         hold on;
%     elseif module_list(modules,5) >0.4
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'b')
%         hold on;
%     end
% end
% 
% %plot scatterplot of module diff (:,3) vs angle diff (:,8) by color depending on same
% %or diff max peak location (:,4)
% 
% 
% figure;
% 
% for modules=1:length(module_list);
%     if module_list(modules,4) == 2 % same max peak wall or corner 
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'r')
%         hold on;
%     elseif module_list(modules,4) == 1 %corner/border same one border
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'g')
%         hold on;
%     elseif module_list(modules,4) == 0 %different wall
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'b')
%         hold on;
%    elseif module_list(modules,4) == 5 %one is centered
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'm')
%         hold on;
%     end
% end
% 
% %plot scatterplot of module diff (:,3) vs angle diff (:,8) by color depending on same
% %or diff max peak location distance (:,5)
% 
% 
% figure;
% 
% for modules=1:length(module_list);
%     if module_list(modules,4) == 2 % same max peak wall or corner 
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'r')
%         hold on;
%     elseif module_list(modules,4) == 1 %corner/border same one border
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'g')
%         hold on;
%     elseif module_list(modules,4) == 0 %different wall
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'b')
%         hold on;
%    elseif module_list(modules,4) == 5 %one is centered
%         scatter(module_list(modules,3), module_list(modules,8), 'markerfacecolor', 'm')
%         hold on;
%     end
% end
% 
% % look at difference in angle depending on same or diff module
% 
% count1=1;
% count2=1;
% 
% for modules=1:length(module_list);
%     if module_list(modules,3) > 0.7
%         same_module_angle(count1)= module_list(modules,8)
%         count1=count1+1;
%     elseif module_list(modules,3) <0.7
%         diff_module_angle(count2)=module_list(modules,8)
%         count2=count2+1;
%     end
% end
% 
% figure; hist(same_module_angle);
% figure; hist(diff_module_angle);
% 
% % look at difference in max peak vs angle depending on same or diff module
% 
% count1=1;
% count2=1;
% 
% for modules=1:length(module_list);
%     if module_list(modules,3) > 0.7 % same module
%         if module_list(modules,8)< 2  % same angle 
%             same_same(count1)= module_list(modules,6)
%             count1=count1+1;
%         elseif module_list(modules,8) > 2 %diff angle
%             diff_same(count2)=module_list(modules,6)
%             count2=count2+1;
%     end
% end
% 
% figure; hist(same_same);
% figure; hist(diff_same);
% 
% 
% 
% disp('')
% 
% % scatter plot of orientation and module with same or diff distance in different colors
%         
%     