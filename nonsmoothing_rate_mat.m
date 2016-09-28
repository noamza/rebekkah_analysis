parms.dir_load_data = 'C:\Users\Dori\Desktop\Rebekkah data\data from non HD PF at least 14 cells'; 

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
%midpoint1 (+pi/2),midpoint2(-pi/2)
parms.num_of_direction_bins=120;
parms.bin_size=3;
parms.sigma = 1;
  
dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;
 
% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};    
    dat = load(file_name); 
    Cell= dat.S;
  
    if ~isempty(Cell.pos.x2)
        dt = Cell.pos.t(2)-Cell.pos.t(1);
        
        % calculate the the rat's head direction (using average of both leds)
        pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t); 
        
        % get rate matrix - using function: Creat_Rate_Map
        non_smooth_rate_mat=CreateRateMap(Cell.pos.x,Cell.pos.y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
    end 

    %% finds peak firing rate
    
non_smooth_rate_mat(isnan(non_smooth_rate_mat)) = 0;

peak_values_len = 0;

peak_values_list = [];

number_zones =[];
number_zones = (unique(dat.S.peak_zone_mat))';
number_zones(find(number_zones==0)) = [];

number_zones_len = length(number_zones);

for cen = 1:number_zones_len 
    find_spk = find(dat.S.peak_zone_mat == cen);
    if sum(find_spk>length(non_smooth_rate_mat))>0
        disp('')
    end
    
    peak_values_len = peak_values_len+1;
    peak_values_list(peak_values_len)= non_smooth_rate_mat(find_spk);
    
end
%     
    
    %% finds mean firing rate of each field 
    
non_smooth_rate_mat(isnan(non_smooth_rate_mat))= 0;
% 
mean_values_len = 0;

number_zones =[];
number_zones = (unique(dat.S.peak_zone_mat))';
number_zones(find(number_zones==0)) = [];

number_zones_len = length(number_zones);

mean_values_list = [];

find_spk = [];
% 
for cen = 1:number_zones_len;
%               
       find_spk = find(dat.S.number_zone_mat == number_zones(cen));
     if sum(find_spk>length(non_smooth_rate_mat))>0
        disp('')  
     end
     
      mean_values_list(cen) = mean(non_smooth_rate_mat(find_spk));   
end    
        
%%changes zone map to mean firing rate values

new_zone_mat= dat.S.number_zone_mat;

for cen=1:number_zones_len;
new_zone_mat(find(new_zone_mat==number_zones(cen))) = mean_values_list(cen); 
end

%%changes zone map to peak firing rate values

peak_zone_mat= dat.S.number_zone_mat;

for cen=1:number_zones_len;
peak_zone_mat(find(peak_zone_mat==number_zones(cen))) = peak_values_list(cen); 
end

% plot mean values
sorted_values = sort(mean_values_list);

vari_over_mean = var(sorted_values)/mean(sorted_values);

peak_rates_list = dat.S.peak_rates;

old_zone_mat = dat.S.zone_mat;

stats= strcat('stats_',file_name);
save(stats, 'i', 'mean_values_list', 'new_zone_mat', 'peak_rates_list', 'old_zone_mat', 'non_smooth_rate_mat', 'peak_values_list', 'peak_zone_mat') 

end



disp('')