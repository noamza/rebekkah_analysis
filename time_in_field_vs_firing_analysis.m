parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=3;
parms.sigma = 1.5;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
     
    % calculate the the rat's head direction (using average of both leds)
        pos_mean_x=(S.pos.x + S.pos.x2)/2;
        pos_mean_y=(S.pos.y + S.pos.y2)/2;
         
   [time_mat]=CreateTimeMap(pos_mean_x,pos_mean_y,S.pos.t,parms);
    
   time_spent_zone=[];
   
   for zone_num= 1:length(S.sorted_means)
       time_spent_zone(zone_num)=sum(time_mat(S.number_zone_mat==zone_num)); %time spent in each field in seconds
   end
   
  % figure; scatter(time_spent_zone, S.peak_rates);
  % figure; imagesc(time_mat);
  % figure; imagesc(S.rate_mat);
   
   max_zone_num= find(S.peak_rates==max(S.peak_rates));
   ratio_max_time_from_mean= time_spent_zone(max_zone_num)/median(time_spent_zone); %the time spent in max firing field over the mean time spent in each field
  
%    figure; scatter(time_spent_zone, S.peak_rates);
   corr_time_zone= corrcoef(time_spent_zone, S.peak_rates);     %correlation coeff of the time spent in the fields and its firing rates
   
   corr_time_zones(i)= corr_time_zone(2);
   ratio_time_mean(i)= ratio_max_time_from_mean;
   
 disp(''); 
   
end


save('time spent vs firing', 'ratio_time_mean', 'corr_time_zones');