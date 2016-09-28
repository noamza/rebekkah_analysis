dbstop if error

parms.dir_load_data = 'N:\users\dori\Matlab_Programs\muscimol\DB_musc_MEC';
parms.dir_save_data = 'N:\users\rebekkah\mus results';
parms.dir_save_images= 'N:\users\rebekkah\mus images';

cd(parms.dir_load_data);

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\DB*.mat'));
parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=3;
parms.sigma = 2.5;
file_names = {dir_list.name};



for i=1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat;
    
    pos_mean_x=[];
    pos_mean_y=[];
    spk_x=[];
    spk_y=[];
    pos_mean_x=(Cell.db.B(1).pos_data.x1 + Cell.db.B(1).pos_data.x2)/2;
    pos_mean_y=(Cell.db.B(1).pos_data.y1 + Cell.db.B(1).pos_data.y2)/2;
    
    
    spk_x=interp1(Cell.db.B(1).pos_data.t,pos_mean_x,Cell.db.B(1).spike_data.ts);
    spk_y=interp1(Cell.db.B(1).pos_data.t,pos_mean_y,Cell.db.B(1).spike_data.ts);
    
    % if first gridness is above 0.3
    
    rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.db.B(1).pos_data.t,spk_x,spk_y,Cell.db.B(1).spike_data.ts,parms);
    rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
    R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
    gridness= GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
    
    [HD_rate_hist,HD_rayleigh_score,HD_rayleigh_angle]=ComputeHeadDirectionality...
        (Cell.db.B(1).pos_data.t,Cell.db.B(1).pos_data.x1,Cell.db.B(1).pos_data.y1,...
        Cell.db.B(1).pos_data.x2,Cell.db.B(1).pos_data.y2,Cell.db.B(1).spike_data.ts,parms);
    
    
    
    if gridness >= 0.3 & HD_rayleigh_score <= 0.25 & length(Cell.db.B)==3
        
        % make sure recovery session is also grid:
        
        pos_mean_x_r=(Cell.db.B(3).pos_data.x1 + Cell.db.B(3).pos_data.x2)/2;
        pos_mean_y_r=(Cell.db.B(3).pos_data.y1 + Cell.db.B(3).pos_data.y2)/2;
        
        
        spk_x_r=interp1(Cell.db.B(3).pos_data.t,pos_mean_x_r,Cell.db.B(3).spike_data.ts);
        spk_y_r=interp1(Cell.db.B(3).pos_data.t,pos_mean_y_r,Cell.db.B(3).spike_data.ts);
        
        % if first gridness is above 0.3
        
        rate_mat_r=CreateRateMap(pos_mean_x_r,pos_mean_y_r,Cell.db.B(3).pos_data.t,spk_x_r,spk_y_r,Cell.db.B(3).spike_data.ts,parms);
        rate_mat_after_auto_pearson_r = Cross_Correlation(rate_mat_r,rate_mat_r);
        
        try 
            R_outer = FindROuter(rate_mat_after_auto_pearson_r,parms);
        gridness_r= GridnessRadius(rate_mat_after_auto_pearson_r,parms,R_outer,i);
        catch ME
            gridness= -5;
        end
        
       % stability= corr2(rate_mat, rate_mat_r);
        
        if gridness_r >= 0.2 %& stability > 0.3;
            
            
            
            
            % find max_peak location and max firing rate for all 7 arenas
            S=[];
            S=get_zones(rate_mat, S);
            
            
            
            %keep max_index_norm, keep max firing rate
            norm_max_indices(1,1)= S.norm_max_index(1);
            norm_max_indices(1,2)= S.norm_max_index(2);
            
            firing_rates(1) = max(max(rate_mat));   %peak firing rates
            gridness_scores(1)= gridness;
            
            n=5;
            m=2;
            fig=figure;
            subplot(n,m,1)
            imagesc(rate_mat); axis equal; axis off; hold on;
           title(sprintf('peak = %0.1f \n mean = %0.1f', max(max(rate_mat)), nanmean2(rate_mat)));
            subplot(n,m,2)
            imagesc(S.zone_mat);
            title(sprintf('gridness=%0.1f', gridness));
            axis equal;axis off;
            
            
            mean_firing(1)= nanmean2(rate_mat);
            
            session_len= length(Cell.db.B(2).pos_data.t);
            divide_sesh= round(session_len/3);
            
            
            a= 1:divide_sesh;
            b= divide_sesh+1:2*divide_sesh;
            c= 2*divide_sesh+1:session_len;
            
            ...........1st
                
        pos_mean_x=[];
        pos_mean_y=[];
        spk_x=[];
        spk_y=[];
        pos_mean_x=(Cell.db.B(2).pos_data.x1(a) + Cell.db.B(2).pos_data.x2(a))/2;
        pos_mean_y=(Cell.db.B(2).pos_data.y1(a) + Cell.db.B(2).pos_data.y2(a))/2;
        
        
        
        spk_ind_1=[];
        spk_ind_2=[];
        spk_ind_3=[];
        
        tmp = (Cell.db.B(2).spike_data.ts-Cell.db.B(2).pos_data.t(divide_sesh));
        tmp(tmp<0)= Inf;
        [spk_ind_1 spk_ind_1] = min(tmp); %index of closest value
        
        
        tmp = (Cell.db.B(2).spike_data.ts-Cell.db.B(2).pos_data.t(divide_sesh*2));
        tmp(tmp<0)= Inf;
        [spk_ind_2 spk_ind_2] = min(tmp); %index of closest value
        
        
        %tmp = (Cell.db.B(2).spike_data.ts-Cell.db.B(2).pos_data.t(session_len));
        %tmp(tmp<0)= Inf;
        %[spk_ind_3 spk_ind_3] = min(tmp); %index of closest value
        spk_ind_3= length(Cell.db.B(2).spike_data.ts);
        
        % build the axis of location when spikes are made
        spk_x_1=interp1(Cell.db.B(2).pos_data.t(a),pos_mean_x,Cell.db.B(2).spike_data.ts(1:spk_ind_1));
        spk_y_1=interp1(Cell.db.B(2).pos_data.t(a),pos_mean_y,Cell.db.B(2).spike_data.ts(1:spk_ind_1));
        
        
        rate_mat=[];
        rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.db.B(2).pos_data.t(a),spk_x_1,spk_y_1,Cell.db.B(2).spike_data.ts(1:spk_ind_1),parms);
        
        rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
        R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
        gridness= GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
        
        
        
        S=[];
        S=get_zones(rate_mat, S);
        
        norm_max_indices(2,1)= S.norm_max_index(1);
        norm_max_indices(2,2)= S.norm_max_index(2);
        firing_rates(2) = max(max(rate_mat));   %peak firing rates
        gridness_scores(2)= gridness;
        
        
        subplot(n,m,3)
        imagesc(rate_mat); axis equal;axis off;
       title(sprintf('peak = %0.1f \n mean = %0.1f', max(max(rate_mat)), nanmean2(rate_mat)));
        subplot(n,m,4)
        imagesc(S.zone_mat);axis equal;axis off;
        title(sprintf('gridness=%0.1f', gridness));
        
        
          mean_firing(2)= nanmean2(rate_mat);
        
        ...........2nd
            pos_mean_x=[];
        pos_mean_y=[];
        pos_mean_x=(Cell.db.B(2).pos_data.x1(b) + Cell.db.B(2).pos_data.x2(b))/2;
        pos_mean_y=(Cell.db.B(2).pos_data.y1(b) + Cell.db.B(2).pos_data.y2(b))/2;
        
        
        % build the axis of location when spikes are made
        spk_x_2=interp1(Cell.db.B(2).pos_data.t(b),pos_mean_x,Cell.db.B(2).spike_data.ts(spk_ind_1+1:spk_ind_2));
        spk_y_2=interp1(Cell.db.B(2).pos_data.t(b),pos_mean_y,Cell.db.B(2).spike_data.ts(spk_ind_1+1:spk_ind_2));
        
        rate_mat=[];
        rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.db.B(2).pos_data.t(b),spk_x_2,spk_y_2,Cell.db.B(2).spike_data.ts(spk_ind_1+1:spk_ind_2),parms);
        
        rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
        
        try
            R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
            gridness= GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
        catch ME
            gridness= -5;
        end
        
        
        try
            S=[];
            S=get_zones(rate_mat, S);
            
            norm_max_indices(3,1)= S.norm_max_index(1);
            norm_max_indices(3,2)= S.norm_max_index(2);
            firing_rates(3) = max(max(rate_mat));   %peak firing rates
            gridness_scores(3)= gridness;
            
            subplot(n,m,6)
            imagesc(S.zone_mat); axis equal;axis off;
            title(sprintf('gridness=%0.1f', gridness));
        catch ME
            norm_max_indices(3,1)= nan;
            norm_max_indices(3,2)= nan;
            firing_rates(3) = max(max(rate_mat));   %peak firing rates
            gridness_scores(3)= nan;
        end
        
        subplot(n,m,5)
        imagesc(rate_mat); axis equal;axis off;
       title(sprintf('peak = %0.1f \n mean = %0.1f', max(max(rate_mat)), nanmean2(rate_mat)));       
        
          mean_firing(3)= nanmean2(rate_mat);
        
        
        if gridness < 0.3 %2nd session must lose gridness
            
            .........3rd
                pos_mean_x=[];
            pos_mean_y=[];
            pos_mean_x=(Cell.db.B(2).pos_data.x1(c) + Cell.db.B(2).pos_data.x2(c))/2;
            pos_mean_y=(Cell.db.B(2).pos_data.y1(c) + Cell.db.B(2).pos_data.y2(c))/2;
            
            
            % build the axis of location when spikes are made
            spk_x_2=interp1(Cell.db.B(2).pos_data.t(c),pos_mean_x,Cell.db.B(2).spike_data.ts(spk_ind_2+1:spk_ind_3));
            spk_y_2=interp1(Cell.db.B(2).pos_data.t(c),pos_mean_y,Cell.db.B(2).spike_data.ts(spk_ind_2+1:spk_ind_3));
            
            rate_mat=[];
            rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.db.B(2).pos_data.t(c),spk_x_2,spk_y_2,Cell.db.B(2).spike_data.ts(spk_ind_2+1:spk_ind_3),parms);
            
            rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
            
            try
                R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
                gridness= GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
            catch ME
                gridness= -5;
            end
            
            
            try
                S=[];
                S=get_zones(rate_mat, S);
                
                norm_max_indices(4,1)= S.norm_max_index(1);
                norm_max_indices(4,2)= S.norm_max_index(2);
                firing_rates(4) = max(max(rate_mat));   %peak firing rates
                gridness_scores(4)= gridness;
                
                subplot(n,m,8)
                imagesc(S.zone_mat);axis equal;axis off;
                title(sprintf('gridness=%0.1f', gridness));
                
            catch ME
                norm_max_indices(4,1)= nan;
                norm_max_indices(4,2)= nan;
                firing_rates(4) = max(max(rate_mat));   %peak firing rates
                gridness_scores(4)= nan;
            end
            
            subplot(n,m,7)
            imagesc(rate_mat);axis equal;axis off;
            title(sprintf('peak = %0.1f \n mean = %0.1f', max(max(rate_mat)), nanmean2(rate_mat)));
            
              mean_firing(4)= nanmean2(rate_mat);
            
            
            if gridness < 0.3 %last session must have lost gridness
                
                
                
                %...........recovery
                
                
                pos_mean_x=[];
                pos_mean_y=[];
                
                pos_mean_x=(Cell.db.B(3).pos_data.x1 + Cell.db.B(3).pos_data.x2)/2;
                pos_mean_y=(Cell.db.B(3).pos_data.y1 + Cell.db.B(3).pos_data.y2)/2;
                
                % build the axis of location when spikes are made
                spk_x_4=interp1(Cell.db.B(3).pos_data.t,pos_mean_x,Cell.db.B(3).spike_data.ts);
                spk_y_4=interp1(Cell.db.B(3).pos_data.t,pos_mean_y,Cell.db.B(3).spike_data.ts);
                
                
                rate_mat=[];
                
                rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.db.B(3).pos_data.t,spk_x_4,spk_y_4,Cell.db.B(3).spike_data.ts,parms);
                
                rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
                
                try
                    R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
                    gridness= GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
                catch ME
                    gridness= -5;
                end
                
                
                
                try
                    S=[];
                    S=get_zones(rate_mat, S);
                    
                    norm_max_indices(4,1)= S.norm_max_index(1);
                    norm_max_indices(4,2)= S.norm_max_index(2);
                    firing_rates(4) = max(max(rate_mat));   %peak firing rates
                    gridness_scores(4)= gridness;
                    
                    subplot(n,m,10);
                    imagesc(S.zone_mat); axis equal; axis off;
                    title(sprintf('gridness=%0.1f', gridness));
                    
                catch ME
                    norm_max_indices(4,1)= nan;
                    norm_max_indices(4,2)= nan;
                    firing_rates(4) = max(max(rate_mat));   %peak firing rates
                    gridness_scores(4)= nan;
                end
                
                
                subplot(n,m,9)
                imagesc(rate_mat); axis equal; axis off;
                     title(sprintf('peak = %0.1f \n mean = %0.1f', max(max(rate_mat)), nanmean2(rate_mat)));
                
                  mean_firing(5)= nanmean2(rate_mat);
                
                %...................
                % see if remains in same location (cluster score)
                
                distance_matrix=nan(length(norm_max_indices),length(norm_max_indices));
                
                for h=1:length(norm_max_indices)-1;
                    for j= h+1:length(norm_max_indices);
                        distance_matrix(h,j)= Distance(norm_max_indices(h,1), norm_max_indices(h,2),norm_max_indices(j,1), norm_max_indices(j,2));
                    end
                end
                
                cluster_matrix= distance_matrix;
                cluster_matrix(distance_matrix>=0.35)= 0;
                cluster_matrix(distance_matrix<0.35)= 1;
                
                cluster_score= sum(sum(cluster_matrix==1))/ sum(sum(~isnan(cluster_matrix)));
                
                %check for both cases where firing decreases and where remains the same
                cd(parms.dir_save_data);
                name= sprintf('%s', file_name);
                save(name, 'norm_max_indices', 'firing_rates', 'cluster_score', 'gridness_scores', 'mean_firing');
                
                cd(parms.dir_save_images);
                saveas(fig, sprintf('image_%s.jpg', file_name));
                cd(parms.dir_load_data);
                
                
                close all
            end
        end
        end
    end
end