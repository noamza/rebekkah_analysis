function [ output_args ] = Untitled6( input_args )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here



% get all files in folder
% create file_list
% create loop 

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\data set for adaptation analysis- 0.4 grid 0.2 hd with circular'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

figure;

for i=1:length(file_names)
    file_name = file_names{i};
    load(file_name);

  %  parms.bin_size=6;
%simulated_rate_map_ns= CreateRateMapNoSmooth(pos_mean_x, pos_mean_y, Cell.pos.t, spk_x, spk_y, spk_t, parms);


%parms.sigma=1.5;
%simulated_rate_map=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,spk_t,parms);

%S_ns= get_zones_nonsmooth(simulated_rate_map_ns, Cell, Cell.PF_radius);

     S_find= get_zones(S.rate_mat, S);
    
     if S_find.number_of_PF > 4
     
    S.peak_rates= S_find.peak_rates;
    S.sorted_means= S_find.sorted_means;
    
    
    figure; imagesc(S.rate_mat)
     figure; imagesc(S_find.zone_mat)
%     max_peak = max(S.peak_rates);
%     second_peak = S.sorted_means(end-1);
%     min_peak = min(S.peak_rates);
    
    %for all
    
    fano_factor(i) = std(S.peak_rates) / mean(S.peak_rates);
    
%     wo_max = S.peak_rates;
%     wo_max(find(wo_max == max_peak)) = [];
%     coeff_over_var_wo_max(i) = std(wo_max)/mean(wo_max);   
% %     
%     wo_2_max = S.peak_rates;
%     wo_2_max(find(wo_2_max == second_peak)) = [];
%     coeff_wo_2_max(i) = std(wo_2_max)/mean(wo_2_max);
%     
%     wo_min= S.peak_rates;
%     wo_min(find(wo_min == min_peak))= [];
%     coeff_wo_min(i) = std(wo_min)/mean(wo_min); 

     else
         fano_factor(i)=nan;
     end 
end

figure; hist(fano_factor);

nanmean(fano_factor)
nanmedian(fano_factor)

% mean_fano = mean(coeff_over_var);
% mean_fano_wo_max= mean(coeff_over_var_wo_max);
% mean_2_max= mean(coeff_wo_2_max);
% mean_min = mean(coeff_wo_min);
%hist(fano_factor_wo_min(:), 0:0.5:4);
%hold on;

%hist(fano_factor(:), 0:0.5:4);

disp('');

% h = findobj(gca, 'Type','patch');
% set(h(1), 'FaceColor','r', 'EdgeColor','k', 'FaceWidth', 0.5)
% set(h(2), 'FaceColor','b', 'EdgeColor','k')
% title('grid cell firing rate variability'); 