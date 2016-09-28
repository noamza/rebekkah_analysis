function AdaptationAnalysis
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


% get all files in folder
% create file_list
% create loop

dbstop if error

parms.dir_load_data = 'C:\Users\Dori\Desktop\best adaptation results';
% parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\adaptation analysis\adaptation analysis results';
%parms.dir_save_images = 'N:\users\rebekkah\final data smoothed\adaptation analysis\analysis\annulus analysis\trajectories nondownscled- greater lesser';

cd(parms.dir_load_data);

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};



for i=1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    % ADD pos mean x shit here from Main.
    
    parms.bin_size=3;
    parms.sigma = 1.5;
    
    if ~isempty(Cell.pos.x2)
        
        dt = Cell.pos.t(2)-Cell.pos.t(1);
        
        % calculate average position of rat
        pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        pos_mean_x_orig= pos_mean_x;
        pos_mean_y_orig=pos_mean_y;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
        %convert it to rate_mat_inds
        bin_size=3;
        
        [size_x size_y]=size(Cell.rate_mat);
        
        pos_t_len= length(Cell.pos.t);
        
        pos_x_inds=[];
        pos_y_inds=[];
        
        if mod(size_x,2) ==0   % Even  number
            pos_x_inds=round(pos_mean_x/bin_size + (size_x/2));
            spk_x_inds= round(spk_x/bin_size + size_x/2);
        elseif mod(size_x,2) == 1    % odd number
            pos_x_inds=round(pos_mean_x/bin_size + (size_x/2) +0.5);
            spk_x_inds= round(spk_x/bin_size + (size_x/2) +0.5);
        end
        
        if mod(size_y,2) ==0   % Even  number
            pos_y_inds= round(pos_mean_y/bin_size + (size_y/2));
            spk_y_inds= round(spk_y/bin_size + size_y/2);
        elseif mod(size_y,2) == 1    % odd number
            pos_y_inds= round(pos_mean_y/bin_size + (size_y/2 +0.5));
            spk_y_inds= round(spk_y/bin_size + size_y/2 +0.5);
        end
        
        
        pos_x_inds(pos_x_inds<=0)=1;
        pos_y_inds(pos_y_inds<=0)=1;
        pos_x_inds(pos_x_inds>max(size_y))= max(size_y);
        pos_y_inds(pos_y_inds>max(size_x))= max(size_x);
        
        
        spk_x_inds(spk_x_inds<=0)=1;
        spk_y_inds(spk_y_inds<=0)=1;
        spk_x_inds(spk_x_inds>max(size_y))= max(size_y);
        spk_y_inds(spk_y_inds>max(size_x))= max(size_x);
        
        parms.sigma=1.5;
        
        %make sure unsmoothed.
        [rate_mat]=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        
        S= Cell;
        S= get_zones(rate_mat, S);
        
% create annulus zone mat
       [size_x, size_y] = size(rate_mat);
%         
        annulus_zone_mat= S.number_zone_mat;

    max_inds=S.max_inds;

        for cen=1:length(max_inds);
            for fig_i =1:size_x
                for j =1:size_y
                    if Distance(fig_i, j, max_inds(cen,1), max_inds(cen,2)) < S.PF_radius * 0.4 %change this depending on how large you want fields to be
                        annulus_zone_mat(fig_i,j)= 0;
                    end
                end
            end
        end

    %S.number_zone_mat= annulus_zone_mat;

   % number_zone_mat= annulus_zone_mat; 

         number_zone_mat=[];
         number_zone_mat= S.number_zone_mat;
        
        spk_x_inds_yy=interp1(Cell.pos.t,pos_x_inds,Cell.spk.t);
        spk_y_inds_yy=interp1(Cell.pos.t,pos_y_inds,Cell.spk.t);
        
        %create vector of expected firing rate at each timestamp (using
        %rate mat
        
        for time=1:length(Cell.pos.t)
            if isnan(pos_x_inds(time))
                firing_rate_t(time)=nan;
            else
                firing_rate_t(time)= rate_mat(pos_y_inds(time), pos_x_inds(time)); %expected firing rate at each timestamp derived using ratemat
            end
        end
        
        
        for time=1:length(Cell.pos.t)
            if isnan(pos_x_inds(time))
                place_field_t(time)=nan;
            else
                place_field_t(time)= number_zone_mat(pos_y_inds(time), pos_x_inds(time)); %number of zonemat at each timestamp
            end
        end
        
        for spk=1:length(spk_x)
            if isnan(spk_x_inds(spk))
                disp('huh')
            else
                spk_zone(spk)= number_zone_mat(spk_y_inds(spk), spk_x_inds(spk)); %
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
        
        
        
        
        spikes_t= spikes_t * dt;   %spikes_t in firing rate converted to firing rate from # of spikes
        
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
        %   speed_value=0;
        %3 second difference
        
        
        %         figure; plot(Cell.pos.t,spikes_t); hold on;
        %          plot(Cell.pos.t,firing_rate_t, 'color', 'r'); hold on;
        %       plot(Cell.pos.t, place_field_t, 'color', 'g');
        
        %  spikes_t=spikes_t_conv;
        
        
        %find instantaeous speed
        dt= 0.02;
        
        for h= 1:length(Cell.pos.t)-1;
            speed= sqrt(((pos_mean_x_orig(h+1)-pos_mean_x_orig(h))^2) + ((pos_mean_y_orig(h+1)-pos_mean_y_orig(h))^2));
            speed=speed * dt;
            speed_t(h)=speed;
        end
        
        speed_t(h+1)=nan;
        
        
        place_field_t(isnan(place_field_t))= 0;
        
        passes=[];
        passes=find(place_field_t~=0);
        j=0;
        
        %%%% CONFIRMED WORKS.
        
        for h=1:length(passes)-2
            diff_of_passes=place_field_t(passes(h))-place_field_t(passes(h+1)); %zero means same field number
            pass_loc= passes(h+1)- passes(h);    %one means same timestamp
            if diff_of_passes==0 && pass_loc==1; %if same field and at same timestamp (aka didnt enter field, leave, then enter again)
                j=j+1;
                pass_predicted_values(count)=pass_predicted_values(count)+ firing_rate_t(passes(h));%+firing_rate_t(passes(h+1));
                pass_actual_values(count)= pass_actual_values(count)+spikes_t(passes(h));%+spikes_t(passes(h+1));
                time_bin_value(count)=j;
                place_field_num(count)= place_field_t(passes(h));
                %   speed_value(count)= speed_value(count) + speed_t(passes(h));
                %   speed_value(count+1)=0;
                pass_predicted_values(count+1)=0;
                pass_actual_values(count+1)=0;
                diff_of_passes2=place_field_t(passes(h+1))-place_field_t(passes(h+2));
                
                end_pass(length(time_bin_value))= passes(h);
                
                %where the pass starts in pos_t intervals
                if j==1
                    start_pass(length(time_bin_value))=passes(h);
                end
                
                if diff_of_passes2~=0        %why is this here? doesnt work without it, but don't remember why
                    
                    pass_predicted_values(count)= pass_predicted_values(count)+firing_rate_t(passes(h+1));
                    pass_actual_values(count)= pass_actual_values(count)+spikes_t(passes(h+1));
                    % speed_value= speed_value(count)+speed_t(passes(h+1));
                    time_bin_value(count)=j+1;
                    
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
    
    
    
    
    
    %  speed_value(end)=[];
    pass_predicted_values(end)= [];
    pass_actual_values(end)=[];
    
    %  speed_value= speed_value./time_bin_value;
    
    predicted_spikes= pass_predicted_values;     %expected NUMBER OF SPIKES
    actual_spikes= pass_actual_values;             %actual NUMBER OF SPIKES
    
    time_bin_value = time_bin_value* 0.02; %value in seconds instead of timestamps
    
    predicted_firing_rate= predicted_spikes./time_bin_value;     %expected firing rate
    actual_firing_rate= actual_spikes./time_bin_value;             %actual firing rate
    
    
    % predicted_firing_rate= pass_predicted_values.*time_bin_value;
    % actual_firing_rate= pass_actual_values.*time_bin_value;
    
    
    % look at overdispersion of max field vs rest of fields
    
    max_zone_num= find(S.peak_rates==(max(S.peak_rates)));
    
    place_field_num_inds= find(place_field_num== max_zone_num);
    
    predicted_minus_actual_rate= predicted_firing_rate- actual_firing_rate;
    
    non_max_field_values= predicted_minus_actual_rate;
    non_max_field_values(place_field_num_inds)= [];
   
    max_field_values= predicted_minus_actual_rate(place_field_num_inds);
%     
    ratio_overdispersion= nanmean(max_field_values)/nanmean(non_max_field_values);
    ratio_overdispersions(i)= ratio_overdispersion;
    median_overdispersion= nanmedian(max_field_values)/nanmedian(non_max_field_values);
    ratio_median_overdispers(i)=median_overdispersion;
    
    
    %  predicted_firing_rate(isnan(predicted_firing_rate))=0;
    
    
    ratio_predicted_actual= actual_firing_rate./predicted_firing_rate;
%     diff_predicted_actual= predicted_firing_rate - actual_firing_rate;
%     sub_divide= (predicted_firing_rate-actual_firing_rate)./predicted_firing_rate;
%     diff_number_spikes= predicted_spikes- actual_spikes;
    
    actual_to_predicted= (actual_firing_rate-predicted_firing_rate)./(actual_firing_rate+predicted_firing_rate);  


%    correlation_ratio= corrcoef(ratio_predicted_actual, time_bin_value,'rows','complete');
%     correlation_difference= corrcoef(diff_predicted_actual, time_bin_value);
%     correlation_sub_divide= corrcoef(sub_divide, time_bin_value);
%     corr_diff_num_spikes= corrcoef(diff_number_spikes, time_bin_value, 'rows','complete');
    
    
correlation_ab= corrcoef(actual_to_predicted, time_bin_value,'rows','complete');

    if length(spk_y) ~= length(spk_y_inds)
        disp('ERROR ERROR')
    end
    
    %       cd(parms.dir_save_data);
    %  save(sprintf('%s',file_name), 'predicted_firing_rate', 'actual_firing_rate', 'time_bin_value',...
    %       'ratio_predicted_actual', 'diff_predicted_actual', 'correlation_ratio', 'correlation_difference',...
    %       'correlation_sub_divide', 'sub_divide', 'diff_number_spikes', 'corr_diff_num_spikes', 'ratio_overdispersion');
    %
    
    
    corr_diffs(i)= correlation_ab(2,1);
    
    
    
    % show trajectories of zeros and non-zeros actual passes
    
    %indices of zeros:
    
   zero_inds= find(actual_firing_rate==0);
    non_zero_inds= find(actual_firing_rate);
    
    greater=find(actual_to_predicted>=0.25);
    lesser= find(actual_to_predicted>=-0.5 & actual_to_predicted <=0);

%greater= non_zero_inds;
%lesser= zero_inds;

    %ADD IF YOU WANT TO DOWNSCALE
    
%     [zero_inds, non_zero_inds] = Downscale(zero_inds, non_zero_inds, Cell.number_of_PF, place_field_num);  
%     
%     traject_nonfire_x = nan(size(pos_mean_x'));
%     traject_nonfire_y = nan(size(pos_mean_x'));
%     
%     traject_fire_x = nan(size(pos_mean_x'));
%     traject_fire_y = nan(size(pos_mean_x'));
%     
    for h=1:length(lesser);
        start= start_pass(lesser(h));
        last= end_pass(lesser(h));
        traject_nonfire_x(start:last)= pos_mean_x(start:last);
        traject_nonfire_y(start:last)= pos_mean_y(start:last);
    end
    
    for h=1:length(greater);
        start= start_pass(greater(h));
        last= end_pass(greater(h));
        traject_fire_x(start:last)= pos_mean_x(start:last);
        traject_fire_y(start:last)= pos_mean_y(start:last);
    end
%     
    %remove not in field values (aka where field_pass_t == 0)
    
    field_zero_inds= find(place_field_t == 0);
    %
    %     traject_non_zeros_x(field_zero_inds) = nan;
    %     traject_non_zeros_y(field_zero_inds) = nan;
    %     traject_zeros_x(field_zero_inds)= nan;
    %     traject_zeros_y(field_zero_inds)= nan;
    
   % doesnt_fire = find(actual_firing_rate==0);
   % sum_zeros=length(doesnt_fire);
    
    
    
    n=4;
    m=2;
  %  cd(parms.dir_save_images);
   % fig = figure;
    %subplot(n,m,[1,2,3,4])
    subplot(n,m,i)
    scatter(time_bin_value, ratio_predicted_actual); 
    %title(sprintf('number of PF= %d, size_of_PF= %f', S.number_of_PF, S.PF_radius));
    xlabel('length of pass (seconds)', 'fontsize', 14, 'fontname', 'calibri');
    ylabel('actual/predicted rate', 'fontsize', 14, 'fontname', 'calibri');
    %axis square;
    %[BestFit,~]=fit(time_bin_value,ratio_predicted_actual,'poly2');
    %plot(BestFit, time_bin_value, ratio_predicted_actual);
    
    % %
    %
    %
    
    %
%     subplot(n,m,5)
%     plot(traject_nonfire_x,traject_nonfire_y, '.',  'MarkerSize', 0.5);hold on;
%     axis ij; axis equal; axis off;
%     title('LESS THAN EXPECTED')
    
%convert values to position coordinates
% sizee= max(pos_mean_x) - min(pos_mean_x);    %only take x since arena is square or circular 
% 
% conv_max_inds= S.max_inds*bin_size - (sizee/2); %convert from rate_mat coords to plot coords
% 
% % 
%     PF_radius_cnvrt= S.PF_radius * 2.7;
%     
%     for h=1:length(S.max_inds)
%         circle(conv_max_inds(h,2),conv_max_inds(h,1), PF_radius_cnvrt)
%     end
%     
%     subplot(n,m, 6)
%     plot(traject_fire_x,traject_fire_y, '.', 'MarkerSize', 0.5);hold on;
%  %   plot(spk_x,spk_y,'.r', 'MarkerSize', 0.5);  %add spikes to plot
%     axis ij;
%     axis equal; axis off;
%     title('GREATER THAN EXPECTED')
%     %
%     
%     for h=1:length(S.max_inds)
%         circle(conv_max_inds(h,2),conv_max_inds(h,1), PF_radius_cnvrt)
%     end
%     
%     %
%     %
%     subplot(n,m,7)
%     imagesc(Cell.zone_mat); hold on; axis equal; axis off;
% %    title(sprintf('didnt fire x= %d total_passes=%d', sum_zeros, length(actual_firing_rate)));
%     subplot(n,m,8);
%     plot(pos_mean_x,pos_mean_y,'k');hold on;
%     plot(spk_x,spk_y,'.r');  %add spikes to plot
%     axis off;  axis ij; axis equal;
%     
%     
%     

%     [med_fired, mean_fired, med_nonfired, mean_nonfired]=  FindClosestPtCOM(start_pass, end_pass, place_field_t, pos_mean_x, pos_mean_y, ...
%                                field_zero_inds, lesser, greater, S.max_inds, S.PF_radius)
%     
%     
%     total_mean_fire(i)= mean_fired;
%     total_med_fire(i)= med_fired;
%     total_mean_nonfire(i)= mean_nonfired;
%     total_med_nonfire(i)= med_nonfired;
%     
%     
%     cd(parms.dir_save_images);
%     saveas(fig, sprintf('%s.jpg', file_name));
%     saveas(fig, sprintf('%s.fig', file_name));
%     %
%     cd(parms.dir_load_data);
    
    
    
   % close all
    
    clearvars -except i file_names parms corr_diffs ratio_overdispersions ratio_median_overdispers total_mean_fire total_med_fire total_mean_nonfire total_med_nonfire
    
    
    
    
end

% difference_mean= total_mean_nonfire- total_mean_fire;
% difference_med= total_med_nonfire- total_med_fire;
% 
% save('trajectory distances from center DOWNSCALE', 'total_mean_fire', 'total_med_fire', 'total_mean_nonfire', 'total_med_nonfire', ...
%                'difference_mean', 'difference_med');






end
