parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
%parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\data sets\images for adaptation analysis nonsmooth';
%parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\data sets\results for adaptation analysis nonsmooth';

parms.bin_size=3;
parms.sigma = 1.5;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

%initialize variables
max_var_z= nan(1, length(file_names));
mean_other_var_z= nan(1, length(file_names));
corrs=nan(1,length(file_names));
corrs_wo_hyper=nan(1,length(file_names));
corrs_hyper=nan(1,length(file_names));

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
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    %can't find how Giocomo paper defined radius, will stick to my 70%
    autocorr=Cross_Correlation(Cell.rate_mat,Cell.rate_mat);
    auto_max_inds = FindAutoMaxInds(autocorr);
    
    [pos_i, pos_j]= ConvertCoordinates(Cell.rate_mat, parms.bin_size,pos_mean_x,pos_mean_y);
    [spk_i, spk_j]= ConvertCoordinates(Cell.rate_mat, parms.bin_size,spk_x,spk_y);
    
    inds_x_y= nan(length(Cell.max_inds),2);
    inds_x_y(:,1)= Cell.max_inds(:,2)* parms.bin_size + min(pos_mean_x);
    inds_x_y(:,2)= Cell.max_inds(:,1)* parms.bin_size + min(pos_mean_y);
    
    max_inds=inds_x_y;
    
    SDM= nan(1,length(spk_x));
    for h=1:length(spk_x)
        dist=nan(1,length(max_inds));
        if ~isnan(spk_x(h))    %if spk isn't nan
            for cen= 1:length(max_inds)
                dist(cen)=Distance(spk_x(h),spk_y(h), max_inds(cen,1), max_inds(cen,2));
            end
            min_ind=find(dist==min(dist));
            min_ind=min_ind(1); %in case two minimums
            SDM(h)=Distance(spk_x(h),spk_y(h), max_inds(min_ind,1), max_inds(min_ind,2));
        end
    end
    
    r= findPlaceFieldRadius(autocorr, auto_max_inds);
    SDM=SDM/r; %divide by field radius so can be compared to all cells
    
    size_x= [min(pos_mean_x) max(pos_mean_x)];
    size_y= [min(pos_mean_y) max(pos_mean_y)];
    dist_boundary=nan(1,length(Cell.pos.t));
    for h=1:length(Cell.pos.t)
        [~,dist_boundary(h)]= findDistPtToBorder(size_x, size_y, [pos_mean_x(h) pos_mean_y(h)]);
    end
    
    % doesn't work, find finds index off by one
    %     inds= 2:length(dist_boundary)-1;
    %     ex= find(dist_boundary(inds) >= 0.1 & (dist_boundary(inds-1) <0.1 ));
    %     en= find(dist_boundary(inds) < 0.1 & (dist_boundary(inds-1) >=0.1 ));
    
    ex=[];
    en=[];
    count=1;
    count2=1;
    for inds= 2:length(dist_boundary)-1;
        if (dist_boundary(inds) >= 0.1 & (dist_boundary(inds-1) <0.1 ));
            ex(count)= inds;
            count=count+1;
        elseif (dist_boundary(inds) < 0.1 & (dist_boundary(inds-1) >=0.1 ));
            en(count2)= inds;
            count2=count2+1;
        end
    end
    
    en(en < ex(1))=[]; %disregard if enters boundary before exiting
    
    lenn=min(length(ex), length(en));
    time_since_boundary=zeros(1,length(dist_boundary));
    count=1;
    for h=1:lenn
        len= en(h)-ex(h);
        time_since_boundary(ex(h):en(h)-1)= 0.02:0.02:0.02*(len);
        %         if len < 500 % only look at trajectories longer than 10 secs (Giocoma used 60 s but for mice)
        %             shorter_10_en(count)= en(h);
        %             shorter_10_ex(count)= ex(h);
        %             count=count+1;
        %         end
    end
    
    timestamps=nan(1,length(Cell.spk.t));
    for h=1:length(Cell.spk.t)
        timestamp= find((Cell.pos.t>=Cell.spk.t(h)));
        if length(timestamp)>= 1    % sometimes spk.t > Cell.pos.t
            timestamps(h)=timestamp(1);
        end
    end
    
    %remove areas where timestamps is nan (sometimes spk.t is larger than
    %pos.t)
    SDM(isnan(timestamps))= []; 
    timestamps(isnan(timestamps))= []; 
        
    %removes spikes that fired inside border (only looking at trajectories
    %after exited boundary area)
    count=1;
    spks_in_border=[];
    for h=1:length(timestamps)
        if time_since_boundary(timestamps(h))==0
            spks_in_border(count)= h;
            count=count+1;
        end
    end
    
    timestamps(spks_in_border)=[];
    SDM(spks_in_border)=[];
    
    figure; scatter(time_since_boundary(timestamps),SDM);
    lsline;
    
    corr_time= corrcoef(time_since_boundary(timestamps),SDM);
    corrs(i)= corr_time(2);
       
    %repeat with removing hyperfield to see if that's the main contributer
    %to the effect seen
    field_t= FindFieldPerTimeBin(Cell.pos.t, pos_i, pos_j, Cell.number_zone_mat);
    
    max_index= find(Cell.peak_rates==max(Cell.peak_rates));
    hyperfield_times= find(field_t == max_index);
    
    count=1;
    remove=[];
    for h=1:length(hyperfield_times)
        for j=1:length(lenn)
        if hyperfield_times(h) >= ex(j) & hyperfield_times(h) <= en(j)
            remove(count)= j;
            count=count+1;
        end
        end
    end
    
    remove=unique(remove);
    
    %embed vector with 0s where hyperfield is involved
    time_since_boundary_edit=time_since_boundary;
    for h=1:length(remove)
        time_since_boundary_edit(ex(h):en(h)-1)= 0;  
    end
    
    
    
    count=1;
    spks_in_hyper=[];
    for h=1:length(timestamps)
        if time_since_boundary_edit(timestamps(h))==0
            spks_in_hyper(count)= h;    %inds of timestamp that include hyperfield-passed trajectories
            count=count+1;
        end
    end
    
    timestamps_wo_hyper=timestamps;
    timestamps_wo_hyper(spks_in_hyper)=[];
    SDM_wo_hyper=SDM;
    SDM_wo_hyper(spks_in_hyper)=[];
    
    corr_time= corrcoef(time_since_boundary(timestamps_wo_hyper),SDM_wo_hyper);
    corrs_wo_hyper(i)= corr_time(2);
    
    % examine only trajectories that passed hyperfield
    timestamps_hyper=timestamps(spks_in_hyper);
    SDM_hyper= SDM(spks_in_hyper);

    corr_time= corrcoef(time_since_boundary(timestamps_hyper),SDM_hyper);
    corrs_hyper(i)= corr_time(2);
    
    %find trajectories that start leaving boundary areato reentering
    % define border area as 0.1 normalized distance
    
    
    %         traject_pos_x=cell(1);
    %         traject_pos_y=cell(1);
    %         traject_ex_firing=cell(1);
    %         count=1;
    %         trajectory_x=[];
    %         trajectory_y=[];
    %         ex_firing=[];
    %         for h=2:length(Cell.pos.t)
    %             if dist_boundary(h) >= 0.1
    %                 trajectory_x(end+1)= pos_mean_x(h);
    %                 trajectory_y(end+1)= pos_mean_y(h);
    %                 ex_firing(end+1)= Cell.rate_mat(pos_i(h), pos_j(h));
    %             elseif dist_boundary(h) < 0.1
    %                 traject_pos_x{count}= trajectory_x;
    %                 traject_pos_y{count}= trajectory_y;
    %                 traject_ex_firing{count}= ex_firing;
    %                 if dist_boundary(h-1)>= 0.1
    %                     count=count+1;
    %                     trajectory_x=[];
    %                     trajectory_y=[];
    %                     ex_firing=[];
    %                 end
    %             end
    %
    %         end
    
    %
    %
    %
    %     beta= 3/2; %beta*r= outer annulus
    %     outer_r= r * beta;
    %
    %     %A=mean peak firing rate over all fields (here we use it as peak firing rate of specific field);
    %     peak_rate= peak_rates(max_ind);
    %     firing_at_x= peak_rate * exp((-x^2)/(2*r^2));
    %
    %
    %     %entry= outer ring to inner, exit= inner ring to outer
    %     %alpha= entry, beta= exit
    %     %angle_phi= angle of entrance and exit trajectory
    %     %x_drift= drift displacement
    %     %lambda_expansion= grid field expansion
    %
    %     prob_en= A
    %
    %
    %
    %     expected_t= FindFiringRatePerTimeBin(Cell.pos.t, pos_i, pos_j, Cell.rate_mat);
    %     dt= Cell.pos.t(2)- Cell.pos.t(1);
    %     expected_t= expected_t*dt; %changes to spikes per timebin instead of per sec
    %
    %     observed_t= FindSpikesPerTimeBin(Cell.pos.t, Cell.spk.t);
    %     %smooth observed:
    %     Win=hamming(13);
    %     Win= Win/sum(Win);
    %     observed_t=conv(observed_t,Win,'same');
    %
    %     PF_radius= findPlaceFieldRadius(Cell.autocorr, Cell.auto_max_inds);
    %     [~, number_zone_mat]= CreateZoneMat(Cell.rate_mat, PF_radius, Cell.max_inds, Cell.peak_rates);
    %
    %     %expected and observed both in NUMBER of spikes (not per sec)
    %     field_t= FindFieldPerTimeBin(Cell.pos.t, pos_x_inds, pos_y_inds, number_zone_mat);
    %
    %     z_t= nan(1,length(Cell.pos.t));
    %     for len= 1:length(Cell.pos.t);
    %         z_t(len)= (observed_t(len)- expected_t(len))/(sqrt(expected_t(len)));
    %     end
    %
    %     %     for h=2:length(field_t);
    %     %         if field_t(h)~=0 & field(h-1)==0
    %     %             trajectory_observed{count}=
    %     %             trajectory_expected{count}=
    %     %
    %     %
    %     %         end
    %     %     end
    %     %
    %
    %     overall_rate= nanmean2(expected_t);
    %
    %     %exp_too_low= find(expected_t <= overall_rate);  %remove intervals where expected_t is too low
    %     %z_t(exp_too_low)= nan;
    %
    %     var_of_z_by_field=nan(1,length(Cell.peak_rates));
    %     for len=1:length(Cell.peak_rates)
    %         var_of_z_by_field(len)= nanvar(z_t(field_t==len));
    %     end
    %
    %
    %     %max_var_z(i)= var_of_z_by_field(max_field_num);
    %     var_of_z_by_field_r= var_of_z_by_field;
    %     %var_of_z_by_field_r(max_field_num)= nan;
    %     %mean_other_var_z(i)= nanmean(var_of_z_by_field_r);
    %
    %     % pivot_score= max(Cell.peak_rates)/ mean(Cell.peak_rates);
    %
    %
    %     %if max(Cell.peak_rates) < 20
    %     %if pivot_score >= 1.8
    %     %max_var_z_pp(count) = var_of_z_by_field(max_field_num);
    %     %var_of_z_by_field_pp_r= var_of_z_by_field;
    %     %var_of_z_by_field_pp_r(max_field_num)= nan;
    %     %mean_other_var_z_pp(count)= nanmean(var_of_z_by_field_pp_r);
    %
    %     n=1;
    %     m=2;
    %     %figure;
    %     %subplot(n,m,2);
    %
    %     %norm_peak_rates= Cell.peak_rates/max(Cell.peak_rates);
    %     %norm_var_of_z_by_field=var_of_z_by_field/min(var_of_z_by_field);
    %
    %     %norm_peak_rates_wo_max= Cell.peak_rates;
    %     %norm_peak_rates_wo_max(max_field_num)= nan;
    %     [size_x size_y]= size(Cell.rate_mat);
    %
    %     %if Cell.max_peak_distance/size_x <= 0.1 & Cell.peak_rates(max_field_num) < 50
    %
    %     %peak_rates_wo_max= Cell.peak_rates;
    %     %peak_rates_wo_max(max_field_num)=nan;
    %
    %     % figure;
    %     % scatter(peak_rates_wo_max, var_of_z_by_field_r); hold on;
    %     % scatter(Cell.peak_rates(max_field_num), var_of_z_by_field(max_field_num), 'FaceColor', 'r')
    %     % title(sprintf('%0.2f', pivot_score));
    %
    %     corrr= corrcoef(Cell.peak_rates, var_of_z_by_field_r);
    %     corrs_all(i)=corrr(2);
    %
    %     % max_var(count)= var_of_z_by_field(max_field_num);
    %     % var_of_z_r= var_of_z_by_field;
    %     % var_of_z_r(max_field_num)= nan;
    %     % mean_other_var_z(count)= nanmean(var_of_z_r);
    %     % subplot(n,m,1); imagesc(Cell.rate_mat);
    %     %end
    %     count=count+1;
    %
    %     %end
    %     %count=count+1;
    %
    %     disp('');
    %
    %
end

figure;
hist(corrs);

mean(corrs)
std(corrs)/sqrt(length(corrs))

mean(corrs_wo_hyper)
std(corrs_wo_hyper)/sqrt(length(corrs_wo_hyper))

mean(corrs_hyper)
std(corrs_hyper)/sqrt(length(corrs_hyper))

% figure; hist(corrs_all);
% mean(corrs_all)

% save('var_z', 'max_var_z', 'mean_other_var_z', 'max_var_z_pp', 'mean_other_var_z_pp');
%Overdispersion in turn is the variance of the z distribution for a set of passes