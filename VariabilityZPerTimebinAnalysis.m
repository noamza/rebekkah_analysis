parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
%parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\data sets\images for adaptation analysis nonsmooth';
%parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\data sets\results for adaptation analysis nonsmooth';

parms.bin_size=3;
parms.sigma = 1.5;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

max_var_z= nan(1, length(file_names));
mean_other_var_z= nan(1, length(file_names));

figure; 
count=1;
% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
      
    % calculate the the rat's head direction (using average of both leds)
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
%     spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
%     spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    

%find z for every time the rat is in field
%compare max field z var to all other fields z var
%or compare z var of all fields

max_field_num= find(Cell.peak_rates==(max(Cell.peak_rates)));

bin_size=3;
[pos_x_inds, pos_y_inds]= ConvertCoordinates(Cell.rate_mat, bin_size,pos_mean_x,pos_mean_y);

expected_t= FindFiringRatePerTimeBin(Cell.pos.t, pos_x_inds, pos_y_inds, Cell.rate_mat);
dt= Cell.pos.t(2)- Cell.pos.t(1);
expected_t= expected_t*dt; %changes to spikes per timebin instead of per sec

observed_t= FindSpikesPerTimeBin(Cell.pos.t, Cell.spk.t);
%smooth observed:
Win=hamming(13);
Win= Win/sum(Win);
observed_t=conv(observed_t,Win,'same');
    
PF_radius= findPlaceFieldRadius(Cell.autocorr, Cell.auto_max_inds); 
[~, number_zone_mat]= CreateZoneMat(Cell.rate_mat, PF_radius, Cell.max_inds, Cell.peak_rates);

%expected and observed both in NUMBER of spikes (not per sec)
field_t= FindFieldPerTimeBin(Cell.pos.t, pos_x_inds, pos_y_inds, number_zone_mat); 

z_t= nan(1,length(Cell.pos.t));
for len= 1:length(Cell.pos.t);
z_t(len)= (observed_t(len)- expected_t(len))/(sqrt(expected_t(len)));   
end


overall_rate= nanmean2(expected_t);

%exp_too_low= find(expected_t <= overall_rate);  %remove intervals where expected_t is too low
%z_t(exp_too_low)= nan;

var_of_z_by_field=nan(1,length(Cell.peak_rates)); 
for len=1:length(Cell.peak_rates)
var_of_z_by_field(len)= nanvar(z_t(field_t==len));
end


%max_var_z(i)= var_of_z_by_field(max_field_num);
var_of_z_by_field_r= var_of_z_by_field;
%var_of_z_by_field_r(max_field_num)= nan;
%mean_other_var_z(i)= nanmean(var_of_z_by_field_r);

% pivot_score= max(Cell.peak_rates)/ mean(Cell.peak_rates);


%if max(Cell.peak_rates) < 20
%if pivot_score >= 1.8
%max_var_z_pp(count) = var_of_z_by_field(max_field_num);
%var_of_z_by_field_pp_r= var_of_z_by_field;
%var_of_z_by_field_pp_r(max_field_num)= nan;
%mean_other_var_z_pp(count)= nanmean(var_of_z_by_field_pp_r);

n=1;
m=2;
%figure; 
%subplot(n,m,2);

%norm_peak_rates= Cell.peak_rates/max(Cell.peak_rates); 
%norm_var_of_z_by_field=var_of_z_by_field/min(var_of_z_by_field);

%norm_peak_rates_wo_max= Cell.peak_rates;
%norm_peak_rates_wo_max(max_field_num)= nan; 
[size_x size_y]= size(Cell.rate_mat);

%if Cell.max_peak_distance/size_x <= 0.1 & Cell.peak_rates(max_field_num) < 50

%peak_rates_wo_max= Cell.peak_rates;
%peak_rates_wo_max(max_field_num)=nan;

% figure;
% scatter(peak_rates_wo_max, var_of_z_by_field_r); hold on;
% scatter(Cell.peak_rates(max_field_num), var_of_z_by_field(max_field_num), 'FaceColor', 'r')
% title(sprintf('%0.2f', pivot_score));

corrr= corrcoef(Cell.peak_rates, var_of_z_by_field_r); 
corrs_all(i)=corrr(2); 

% max_var(count)= var_of_z_by_field(max_field_num);
% var_of_z_r= var_of_z_by_field;
% var_of_z_r(max_field_num)= nan;
% mean_other_var_z(count)= nanmean(var_of_z_r);
% subplot(n,m,1); imagesc(Cell.rate_mat); 
%end
count=count+1;

%end
%count=count+1;

disp('');


end

figure; hist(corrs_all);
mean(corrs_all)

% save('var_z', 'max_var_z', 'mean_other_var_z', 'max_var_z_pp', 'mean_other_var_z_pp');
%Overdispersion in turn is the variance of the z distribution for a set of passes