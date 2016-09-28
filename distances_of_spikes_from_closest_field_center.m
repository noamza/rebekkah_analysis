cd('N:\users\rebekkah\results and info of analysis')
load('variability and border distances.mat')

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

parms.bin_size=3;
parms.sigma=1.5;

spk_dist_hyperfield= nan(1,length(file_names));
spk_dist_nonmax_fields= nan(1,length(file_names));
spk_dist_border= nan(1,length(file_names));
spk_dist_nonborder= nan(1,length(file_names));

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    % calculate the the rat's head direction (using average of both leds)
    pos_mean_x=(S.pos.x + S.pos.x2)/2;
    pos_mean_y=(S.pos.y + S.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(S.pos.t,pos_mean_x,S.spk.t);
    spk_y=interp1(S.pos.t,pos_mean_y,S.spk.t);
    
    rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,S.pos.t,spk_x,spk_y,S.spk.t,parms);
    
    %convert max_inds to spike_mat coordinates
    [pos_x_inds, pos_y_inds] = ConvertCoordinates(rate_mat, parms.bin_size, pos_mean_x,pos_mean_y);
    [spkx_inds, spky_inds] = ConvertCoordinates(rate_mat, parms.bin_size, spk_x, spk_y);
    
    max_inds= max_indices_sm{i};
    distances_spks_from_peak=cell(1,length(max_inds));
    for h=1:length(spk_x)
        if ~isnan(spkx_inds(h))
            distances_from_peaks=nan(1,length(max_inds));
            for cen= 1:length(max_inds)
                distances_from_peaks(cen)= Distance(spkx_inds(h), spky_inds(h), max_inds(cen,1), max_inds(cen,2)); %find distance of spk from zones
            end
            min_dist_zone_num= find(distances_from_peaks==min(distances_from_peaks)); %finds which zone spk belongs to
            min_dist_zone_num=min_dist_zone_num(1); %if two min distances, chooses one
            spk_dist_COM= distances_spks_from_peak{min_dist_zone_num};
            current_len= length(distances_spks_from_peak{min_dist_zone_num});
            spk_dist_COM(current_len+1)=distances_from_peaks(min_dist_zone_num);
            distances_spks_from_peak{min_dist_zone_num}= spk_dist_COM;
        end
    end
    
    peak_rates=peak_rates_all_sm{i};
    mean_spk_dist_zone= nan(1,length(peak_rates));
    for cen=1:length(peak_rates)
        spk_dists= distances_spks_from_peak{cen};
        mean_spk_dist_zone(cen)= mean(spk_dists);
    end
    
    [~, ~, border_zones] = BorderNonborderDiff(rate_mat, max_inds, peak_rates);
    
    mean_spk_dist_border= mean_spk_dist_zone(border_zones);
    mean_spk_dist_nonborder= setdiff(mean_spk_dist_zone, mean_spk_dist_border);
    
    
    
    % max_ind= find(peak_rates==max(peak_rates));
    % mean_spk_dist_zone_wo_max= mean_spk_dist_zone;
    % mean_spk_dist_zone_wo_max(max_ind)=[];
    
    %finds field radius
    autocorr= Cross_Correlation(rate_mat, rate_mat);
    auto_max_inds= FindAutoMaxInds(autocorr);
    PF_radius= findPlaceFieldRadius(autocorr, auto_max_inds);
    
    spk_dist_border(i)= mean(mean_spk_dist_border)/ PF_radius;
    spk_dist_nonborder(i)= mean(mean_spk_dist_nonborder)/ PF_radius;
    
    % divides by PF radius so distances can be comapred between cells
    %     spk_dist_hyperfield(i)= mean_spk_dist_zone(max_ind)/PF_radius;
    %     spk_dist_nonmax_fields(i)= mean(mean_spk_dist_zone_wo_max)/PF_radius;
    
    
    disp('')
end

%found no difference between hyperfield and mean of rest of fields:
%spks do not fire closer to center of hyperfield than other fields.

mean(spk_dist_border)
std(spk_dist_border)/sqrt(86)
mean(spk_dist_nonborder)
std(spk_dist_nonborder)/sqrt(86)

% mean(spk_dist_hyperfield)
% std(spk_dist_hyperfield)/sqrt(86)
% mean(spk_dist_nonmax_fields)
% std(spk_dist_nonmax_fields)/sqrt(86)