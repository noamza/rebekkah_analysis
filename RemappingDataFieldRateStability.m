parms.dir_load_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\max';
parms.dir_save_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\cluster criteria 0.45';
%parms.dir_save_data2 = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results 2 above 0.3';
%parms.dir_save_pictures='N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\images above 0.3';
%parms.dir_save_pictures2='N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\images 2 above 0.3';

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=5;
parms.sigma = 3;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    % finds properties of original arena:
    [size_arena_orig_x, size_arena_orig_y]= size(zone_mats.arena1);
        
    max_inds = FindMaxIndsRateMap(rate_mats.arena1);
    autocorr= Cross_Correlation(rate_mats.arena1, rate_mats.arena1);
    auto_max_inds= FindAutoMaxInds(autocorr);
    PF_radius= findPlaceFieldRadius(autocorr, auto_max_inds);
    
    max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mats.arena1);
    
    %     for h=1:length(max_inds)
    %         if zone_mats.arena1(max_inds(h,1), max_inds(h,2))==0
    %             max_inds(h,1) = nan;
    %             max_inds(h,2) = nan;
    %         end
    %     end
    
    %   max_inds(isnan(max_inds))= [];
    
  [max_inds_len, ~]= size(max_inds);
    
    for  h=1:max_inds_len;
        firing_rates.arena1(h)= zone_mats.arena1(max_inds(h,1), max_inds(h,2));
    end
    
    sorted_means_orig= sort(firing_rates.arena1);
    ordered_inds=nan(1,max_inds_len);
    
            for h=1:length(firing_rates.arena1)
                ordered_inds(h)= find(sorted_means_orig(h) == firing_rates.arena1);
            end
    
    for arena_count= 2:5    % =2:4 if rescaled arena, =5 if last identical arena 
        
        next_arena= zone_mats.(sprintf('arena%d', arena_count));
        
        [size_arena_next_x, size_arena_next_y]= size(next_arena);
        
        %full stretch:
        factor_x= size_arena_orig_x/ size_arena_next_x;
        factor_y= size_arena_orig_y/ size_arena_next_y;
        tform = maketform('affine',[factor_y 0 0; 0 factor_x 0; 0 0 1]);
        %form= [size_square(1), stretch_axis];
        stretched_arena.(sprintf('stretch%d', arena_count)) = imtransform(next_arena,tform);
        
        current_arena=stretched_arena.(sprintf('stretch%d', arena_count));
        
        rates=nan(1,max_inds_len);
        for h=1:max_inds_len
            rates(h)= current_arena(max_inds(h,1),max_inds(h,2));
        end
        
        name=sprintf('arena%d', arena_count);
        
        rates= rates(ordered_inds);
        
        firing_rates.(sprintf('%s', name))= rates;
        
        
    end
    
    figure;
    plot(firing_rates.arena1, 'o-'); hold on;
    plot(firing_rates.arena2, 'o-'); hold on;
    plot(firing_rates.arena3, 'o-'); hold on;
    plot(firing_rates.arena4, 'o-'); hold on;
    plot(firing_rates.arena5, 'o-'); hold on;
    
    
end
