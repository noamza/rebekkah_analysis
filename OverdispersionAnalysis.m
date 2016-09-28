function AdaptationAnalysis

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results for adaptation analysis';
parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\adaptation analysis\results normalized';
parms.dir_save_images = 'N:\users\rebekkah\final data smoothed\adaptation analysis\images normalized';

cd(parms.dir_load_data);

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

parms.bin_size=3;
parms.sigma = 1.5;

for i=1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
       
    dt = Cell.pos.t(2)-Cell.pos.t(1);
    
    % calculate average position of rat
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    [rate_mat]=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
    %convert it to rate_mat_inds
    [pos_x, pos_y]= ConvertCoordinates(rate_mat, parms.bin_size, pos_x,pos_y);
    [spk_x_inds, spk_y_inds]= ConvertCoordinates(rate_mat, parms.bin_size, pos_x,pos_y);
    
    %Create autocorr and find max inds
    autocorr= Crodd_Correlation(rate_mat, reat_mat);
    auto_max_inds= FindAutoMaxInds(autocorr);
    PF_radius= findPlaceFieldRadius(autocorr, auto_max_inds);
    max_inds=FindMaxIndsRateMap(rate_mat);
    
    %find peak rates of each field
    peak_rates= nan(1,length(max_inds));
    for cen= 1:length(max_inds);
        peak_rates(cen)= rate_mat(max_inds(cen,1), max_inds(cen,2));
    end
    
    %create number zone mat
    [~, number_zone_mat]=CreateZoneMat(rate_mat, PF_radius, max_inds, peak_rates);
    
    %create vector of expected firing rate and field location at each timestamp (using
    %rate mat
    
    firing_rate_t= nan(1,length(Cell.pos.t));
    place_field_t= nan(1,length(Cell.pos.t));
    for time=1:length(Cell.pos.t)
        if isnan(pos_x_inds(time))
            firing_rate_t(time)=nan;
            place_field_t(time)=nan;
        else
            firing_rate_t(time)= rate_mat(pos_x(time), pos_y(time)); %expected firing rate at each timestamp derived using ratemat
            place_field_t(time)= number_zone_mat(pos_x(time), pos_y_inds(time));
        end
    end
    
    
    for spk=1:length(spk_x)
        if isnan(spk_x_inds(spk))
            disp('huh')
        else
            spk_zone(spk)= number_zone_mat(spk_x_inds(spk), spk_y_inds(spk)); %
        end
    end
    
    %%firing_rate_t is expecting firing rate
    %%spikes_t is actual firign rate
    
    
    
    %whenver during time place_field enters new field and exits new
    %field, add number of expected firing (from firing rate t) and
    %actual firing rate (from Cell.spk.t and pos.t)
    spikes_t=zeros(1, length(Cell.pos.t));
    %%
    %        train=zeros(1,max(Cell.pos.t));
    %        spk_round=round((Cell.spk.t)*1000);
    %        train(spk_round)=1;
    %        Win=hamming(3);
    %        train=conv(train,Win,'same');
    
    
    
    %
    for h= 1:length(Cell.spk.t)
        tmp=abs(Cell.pos.t-Cell.spk.t(h));
        [~, idx] = min(tmp);
        closest=Cell.pos.t(idx);
        spikes_t(Cell.pos.t==closest)= spikes_t(Cell.pos.t==closest)+1;
    end
    
    
    
    
    %spikes_t= spikes_t/dt;   %spikes_t in firing rate converted to firing rate from # of spikes
    
    %       spikes_t=smooth(spikes_t, 23);
    %       firing_rate_t=smooth(firing_rate_t);
    
    %spikes_t=spikes_t/dt;
    
    %   Win=hamming(11);
    %   spikes_t=conv(spikes_t,Win,'same');
    % firing_rate_t= conv(firing_rate_t, Win,'same');
    
    
    %  spikes_t=spikes_t/dt;
    %  spikes_t_smooth=spikes_t_smooth/dt;
    % spikes_t_conv=spikes_t_conv/dt;
    
    firing_rate_t= firing_rate_t * dt; %number of spikes per timestamp instead of firing rate
    
    count=1;
    pass_predicted_values=0;
    pass_actual_values=0;
    time_bin_value=0;
    
    %3 second difference
    
    
    %         figure; plot(Cell.pos.t,spikes_t); hold on;
    %          plot(Cell.pos.t,firing_rate_t, 'color', 'r'); hold on;
    %       plot(Cell.pos.t, place_field_t, 'color', 'g');
    
    %  spikes_t=spikes_t_conv;
    
    passes=find(place_field_t~=0);
    j=0;
    for h=1:length(passes)-2
        diff_of_passes=place_field_t(passes(h))-place_field_t(passes(h+1));
        if diff_of_passes==0;
            j=j+1;
            pass_predicted_values(count)=pass_predicted_values(count)+ firing_rate_t(passes(h));%+firing_rate_t(passes(h+1));
            pass_actual_values(count)= pass_actual_values(count)+spikes_t(passes(h));%+spikes_t(passes(h+1));
            time_bin_value(count)=j;
            pass_predicted_values(count+1)=0;
            pass_actual_values(count+1)=0;
            diff_of_passes2=place_field_t(passes(h+1))-place_field_t(passes(h+2));
            if diff_of_passes2~=0
                count=count+1;
                j=0;
            end
        else
            % count= count+1; % is this necessary?
            j=0;
        end
        
    end
end

count=1;


%      for h= 1:length(Cell.spk.t);
%         end_of_time_bin = Cell.spk.t(h) + 3;
%         end_of_time_bin_ind= find(Cell.spk.t>end_of_time_bin);
%         end_of_time_bin_ind=end_of_time_bin_ind(1);
%         time_bin_value(count)= Cell.spk.t(end_of_time_bin_ind-1)- Cell.spk.t(h);
%         spikes(count)= (end_of_time_bin_ind-1)- h;
%     %    predicted_rate= rate_mat(pos_mean_x(h), pos_mean_y(h)))
%         h= time_bin_value_ind;
%      end
%






pass_predicted_values(end)= [];
pass_actual_values(end)=[];


predicted_spikes= pass_predicted_values;     %expected NUMBER OF SPIKES
actual_spikes= pass_actual_values;             %actual NUMBER OF SPIKES

time_bin_value = time_bin_value* 0.02; %value in seconds instead of timestamps

predicted_firing_rate= predicted_spikes./time_bin_value;     %expected firing rate
actual_firing_rate= actual_spikes./time_bin_value;             %actual firing rate


% predicted_firing_rate= pass_predicted_values.*time_bin_value;
% actual_firing_rate= pass_actual_values.*time_bin_value;

predicted_firing_rate(isnan(predicted_firing_rate))=0;

ratio_predicted_actual= predicted_spikes./actual_spikes;
diff_predicted_actual= predicted_firing_rate - actual_firing_rate;
sub_divide= (predicted_firing_rate-actual_firing_rate)./predicted_firing_rate;
diff_number_spikes= predicted_spikes- actual_spikes;


correlation_ratio= corrcoef(ratio_predicted_actual, time_bin_value,'rows','complete');
correlation_difference= corrcoef(diff_predicted_actual, time_bin_value);
correlation_sub_divide= corrcoef(sub_divide, time_bin_value);
corr_diff_num_spikes= corrcoef(diff_number_spikes, time_bin_value, 'rows','complete');


if length(spk_y) ~= length(spk_y_inds)
    disp('ERROR ERROR')
end

cd(parms.dir_save_data);
save(sprintf('%s',file_name), 'predicted_firing_rate', 'actual_firing_rate', 'time_bin_value',...
    'ratio_predicted_actual', 'diff_predicted_actual', 'correlation_ratio', 'correlation_difference',...
    'correlation_sub_divide', 'sub_divide', 'diff_number_spikes', 'corr_diff_num_spikes');

cd(parms.dir_save_images);
fig = figure;
scatter(time_bin_value, ratio_predicted_actual);

saveas(fig, sprintf('%s.jpg', file_name));

cd(parms.dir_load_data);

corr_diffs(i)= correlation_ratio(2,1);

close all

clearvars -except i file_names parms corr_diffs

end

figure; hist(corr_diffs);

end
