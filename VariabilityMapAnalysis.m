%% Write function that finds the firing rate variability at each spatial bin

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
%parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\data sets\images for adaptation analysis nonsmooth';
%parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\data sets\results for adaptation analysis nonsmooth';

parms.bin_size=3;
parms.sigma = 1.5;

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
    
    pos_mean_x=[];
    pos_mean_y=[];
    
    % calculate the the rat's head direction (using average of both leds)
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    

    
    % get rate matrix - using function: Creat_Rate_Map
    rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    
    %convert pos_mean_x and pos_mean_y to rate map coords
    bin_size=3;
    [pos_x_inds, pos_y_inds]= ConvertCoordinates(rate_mat, Cell.pos.t,bin_size,pos_mean_x,pos_mean_y);
    
    firing_rate_t=[];
   
    
    spk_idx=interp1(Cell.pos.t,[1:length(Cell.pos.t)],Cell.spk.t);
   
    pos_rate=hist(spk_idx,1:length(Cell.pos.t));
    %should be same as observed_t
    
    
    %FENTON METHOD: 
    %the entire session was divided into 5 s intervals. 
    %For each interval we calculated the expected number of spikes, exp, as:
    %exp= sum(
        
    % create expected firing rate vector at each time bin
    expected_t = FindFiringRatePerTimeBin(Cell.pos.t, pos_x_inds, pos_y_inds, rate_mat); %spikes per sec
       
    %For each selected 5 s interval, 
    %we then calculated z, the normalized standard deviation of obs, 
    %the observed number of spikes as:
    %(obs-exp) / sqrt(exp)
    
    % actual firing rate vector at each time bin
    observed_t=zeros(1, length(Cell.pos.t));    
    
    for h= 1:length(Cell.spk.t) 
        tmp=abs(Cell.pos.t-Cell.spk.t(h));
        [~, idx] = min(tmp);
        observed_t(idx)= observed_t(idx)+1; %in number of spikes per timebin (dt=0.02)
    end
    
     %smooth spikes_t


       


      dt=(Cell.pos.t(2)-Cell.pos.t(1)); 
      observed_t= observed_t / dt; %spikes/sec (instead of spikes/dt)
      pos_rate= pos_rate /dt; %make sure its the same
      
        Win=hamming(13);
        observed_t=conv(observed_t,Win,'same');
        pos_rate= conv(pos_rate, Win,'same');
        
    firing_var_t_check=[];
    var_mat=[];
    
    
    
    
    % vector of actual minus expected values:
    obs_minus_exp=abs(observed_t-expected_t); %% is absolute value necessary
    
    
    
    
    %calculate variability from firing mean at each time point
    %BELOW JUST TO DOUBLE CHECK ABOVE CORRECT
    % for time=1:length(firing_rate_t)
    %     if ~isnan(pos_x_inds(time))
    %    firing_var_t_check(time)=spikes_t(time)- rate_mat(pos_y_inds(time), pos_x_inds(time));
    %     else firing_var_t_check(time)=nan;
    %     end
    % end
    
    
    
    
    % create map of mean variability at each spatial bin
    var_map=CreateVarMap(pos_mean_x,pos_mean_y,Cell.pos.t,parms,obs_minus_exp, rate_mat);
    
 %  figure; imagesc(var_map);
    
    overall_rate= mean2(rate_mat) * 1.1;
    stand_dev= std2(rate_mat);
    
    too_low= stand_dev;
    
  %  remove_coords=find(rate_mat<=too_low);
  %  var_map(remove_coords) =0;
 
 
    
%     n=1; m=2;
%     figure;
%     subplot(n,m,1)
%     imagesc(rate_mat);
%     subplot(n,m,2)
%     imagesc(var_map);
    
  [zone_var_map] =CreateZoneMat(var_map, rate_mat, Cell.max_inds);
    
%     n=1; m=2;
%     figure;
%     subplot(n,m,1)
%     imagesc(Cell.zone_mat);
%     subplot(n,m,2)
%     imagesc(zone_var_map);
    
    
    %find variability at center of fields
    %check fano factors of variability
    var_z_per_field= nan(1,length(Cell.peak_rates));
    for h=1:length(Cell.peak_rates)
        mean_var_of_field= mean(var_map(Cell.number_zone_mat==h));
        var_z_per_field(h)= mean_var_of_field;
        %var_z_per_field_2(h)= var_map(Cell.peak_zone_mat==h);
    end
        
    var_z_fano(i)= var(var_z_per_field)/ mean(var_z_per_field); 
    
    disp('');
    
end

figure; hist(var_z_fano);

