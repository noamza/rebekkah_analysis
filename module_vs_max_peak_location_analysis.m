load('C:\Users\Dori\Desktop\Rebekkah data\module_list_distances_PF_8_halfsmooth')

count=1;
count2=1;
count3=1;
count4=1;

%% add 6th column same or not same based on distance (0.6)
for h= 1:length(module_list);
    if module_list(h,5) <= 0.5 
        module_list (h, 6) = 2;
    elseif module_list(h,5) <0.7 & module_list(h,5) >0.5
        module_list(h,6) = 0;
    elseif module_list(h,5) > 0.7
        module_list (h,6) = 0;
    end
end 

same_module_locations = [];
non_same_module_locations=[];
semi_same_module_locations = [];

for k = 1:length(module_list)
    if module_list(k,3) < 0.8
        non_same_module_locations(count) = module_list(k,6);
        count = count+1
%     elseif module_list(k,3) > 0.6 && module_list(k,3) < 0.8 
%         semi_semi_same_module_locations(count4) = module_list(k,6);
        count4=count4+1;
    elseif module_list(k,3) > 0.8 && module_list(k,3) < 0.9 
        semi_same_module_locations(count2) = module_list(k,6);
        count2=count2+1;
    elseif module_list(k,3) > 0.9
        same_module_locations(count3) = module_list(k,6);
        count3=count3+1;

    end
end 



% same_module_locations(find(same_module_locations==1))= 0;
% semi_same_module_locations(find(semi_same_module_locations==1))= 0;
% % semi_semi_same_module_locations(find(semi_semi_same_module_locations==1))= 0;
% non_same_module_locations(find(non_same_module_locations==1))= 0;


figure;

n= 1
m= 3

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

% subplot(1,4,3)
% bar(hist(semi_semi_same_module_locations) ./ sum(hist(semi_semi_same_module_locations)))
% title('Semi Semi Same')
% ylim([0 0.8])


disp('');
% subplot(1,4,3)
% hist(semi_not_same_module_locations);
% title('Semi Not Same')