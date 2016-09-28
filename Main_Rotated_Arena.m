function Main
% Date: 24 of July 2014
dbstop if error

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\data set for adaptation analysis- 0.4 grid 0.2 hd with circular';
%parms.dir_load_data2 = 'N:\users\rebekkah\final data smoothed\data sets\fitted arena with HD without circular rooms_ 202 cells';
parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\data sets\images for adaptation analysis nonsmooth';
% parms.dir_save_pictures2= 'C:\Users\Dori\Desktop\Rebekkah data\results testing';
% parms.dir_save_pictures3= 'C:\Users\Dori\Desktop\Rebekkah data\results testing';
parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\data sets\results for adaptation analysis nonsmooth';

parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
%midpoint1 (+pi/2),midpoint2(-pi/2)
parms.num_of_direction_bins=120;
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
    %    clear dat;
    %save speed_thresh_hold vector
    % Cell.speed_thresh_hold = speed_thresh_hold;
    
    
    %make sure that new max_inds is used and not old saved one from inaccurate
    %rate_mat and autocorr
    S.max_inds= [];
    S.auto_max_inds=[];
    too_close=[];
    PF_areas=[];
    rate_over_area=[];
    
    
    %% ADD ONLY FOR ROTATING ARENA
    max_pts_inside= find(Cell1.phi_of_arena.percent_inside_final(:,1)==max(Cell1.phi_of_arena.percent_inside_final(:,1)));
    arena_phi= Cell1.phi_of_arena.percent_inside_final(max_pts_inside,5);
    arena_phi= -arena_phi(1);
    
    arena_phi= degtorad(arena_phi); %change to radians
    
    
    %arena_phi=degtorad(0); %for checking
    
    %rotate positions according to arena angle shift
    [pos_x, pos_y]= rotatePath(Cell.pos.x, Cell.pos.y, arena_phi);
    [pos_x2, pos_y2]=rotatePath(Cell.pos.x2, Cell.pos.y2, arena_phi);
    
    Cell.pos.x= pos_x;
    Cell.pos.y= pos_y;
    Cell.pos.x2= pos_x2;
    Cell.pos.y2= pos_y2;
    
    
    
    [spk_x, spk_y]= rotatePath(Cell.spk.x, Cell.spk.y, arena_phi);
    [spk_x2, spk_y2]=rotatePath(Cell.spk.x2, Cell.spk.y2, arena_phi);
    
    Cell.spk.x= spk_x;
    Cell.spk.y= spk_y;
    Cell.spk.x2= spk_x2;
    Cell.spk.y2= spk_y2;
    
    
    if ~isempty(Cell.pos.x2)
        
        dt = Cell.pos.t(2)-Cell.pos.t(1);
        
        % calculate the the rat's head direction (using average of both leds)
        pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
        % get rate matrix - using function: Creat_Rate_Map
        %  rate_mat=CreateRateMap(Cell.pos.x,Cell.pos.y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        [rate_mat]=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        Cell.rate_mat = rate_mat;
        
        Cell.spike_mat = spike_mat;
        
        S=get_zones(rate_mat,dat.S);
        
        % use get matrix after auto pearson correlation - using function: Auto_Pearson_Correlation
        rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
        
        Cell.rate_mat_after_auto_pearson = rate_mat_after_auto_pearson;
        
        R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
        % TODO: to print value of gridness2
        % TODO: to check if the figure print is really required
        [gridness3] = GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
        
        Cell.gridness3 = gridness3;
        pop_gridness2(count)=gridness3;
        
        % calculate hist of HD
        Hist_Moving_Directionality_In_Head_Directionality = ComputePrecentOfMovingDirectionalityInHeadDirectionality(Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2);
        ax = deg2rad(-180:180);
        %hist(Hist_Moving_Directionality_In_Head_Directionality,ax);
        [~, location] = max(hist(Hist_Moving_Directionality_In_Head_Directionality,...
            ax));
        alpha = ax(location);
        Cell.alpha = alpha;
        
        % calculate HD
        [HD_rate_hist,HD_rayleigh_score,HD_rayleigh_angle]=ComputeHeadDirectionality...
            (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2,Cell.spk.t,parms);
        pop_rayleighHD_score(count)=HD_rayleigh_score;
        pop_rayleighHD_angle(count)=HD_rayleigh_angle;
        Cell.HD.rate_hist = HD_rate_hist;
        Cell.HD.rayleigh.score=HD_rayleigh_score;
        Cell.HD.rayleigh.angle=HD_rayleigh_angle;
        
        % calculate the location of the rat by both leds
        pos_average_x = (Cell.pos.x + Cell.pos.x2)/2;
        pos_average_y = (Cell.pos.y + Cell.pos.y2)/2;
        
        % calculate MD
        [MD_rate_hist,MD_rayleigh_score,MD_rayleigh_angle]=ComputeMovingDirectionality...
            (Cell.pos.t,pos_average_x,pos_average_y,Cell.spk.t,parms);
        pop_rayleighMD_score(count) = MD_rayleigh_score;
        pop_rayleighMD_angle(count) = MD_rayleigh_angle;
        Cell.MD.rate_hist = MD_rate_hist;
        Cell.MD.rayleigh.score = MD_rayleigh_score;
        Cell.MD.rayleigh.angle = MD_rayleigh_angle;
        
        % insert variance and mean of gc size
        Cell.variance = S.variance;
        Cell.mean = S.mean_values_list;
        Cell.overall_mean= S.overall_mean;
        Cell.number_of_PF=S.number_of_PF;
        Cell.PF_radius=S.PF_radius;
        Cell.PF_coverage=S.PF_coverage;
        Cell.sorted_means=S.sorted_means;
        Cell.norm_means = S.norm_means;
        Cell.where= S.where;
        Cell.peak_rates = S.peak_rates_list;
        Cell.max_index = S.max_index ;
        Cell.norm_max_index= S.norm_max_index;
        Cell.zone_mat = S.zone_mat;
        Cell.number_zone_mat= S.number_zone_mat;
        Cell.peak_zone_mat= S.peak_zone_mat;
        Cell.max_peak_distance= S.max_peak_distance;
        Cell.max_inds = S.max_inds;
        Cell.location= S.location;
        Cell.i = S.i;
        Cell.auto_max_inds= S.auto_max_inds;
        Cell.rate_mat_after_auto_pearson = S.rate_mat_after_auto_pearson;
        Cell.zone_spike_mat=S.zone_spike_mat;
        Cell.border_to_nonborder_diff_wo_max=S.border_to_nonborder_diff_wo_max;
        Cell.border_to_nonborder_diff=S.border_rate_diff;
        
        Cell_num(count)=i;
        
        % plots and graphs figures
        
        
        bPrintPicture = true;
        if (bPrintPicture)
            % plot all the picturs and graphs
            strTitle = sprintf('Cell_i=%d_r%s_d%s_s%s_t%d_c%d',i,Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell);
            fig =figure('name',strTitle);
            % title(strTitle);
            
            % plot the rat's path (and the spikes on it)
            n = 2;
            m = 5;
            
            subplot(n,m,1);
            plot(pos_mean_x,pos_mean_y,'k');hold on;
            plot(spk_x,spk_y,'.r');
            axis equal;axis off;
            axis ij
            title(file_name);
            
            hold on;
            
            
            
            subplot(n,m,2);
            imagesc(rate_mat); axis image
            
            subplot(n,m,3);
            imagesc(rate_mat_after_auto_pearson); axis image
            hold on;
            title(sprintf('gridness=%0.2f /n elis_grid=%0.2f', S.gridness2, gridness3));
            
            subplot(n,m,4);
            hist(Hist_Moving_Directionality_In_Head_Directionality,ax);
            title('Moving Directionality In Head Directionality');
            
            subplot(n,m,4);
            set(0,'defaulttextinterpreter','none');
            step_len = 360/parms.num_of_direction_bins;
            y_presentation_axis=-180:step_len:180;
            bar(y_presentation_axis,Cell.HD.rate_hist);
            xlabel('angle');
            ylabel('count');
            strTitle = sprintf('Head_Directionality ray=%0.3f',Cell.HD.rayleigh.score);
            title(strTitle);
            
            subplot(n,m,6);
            set(0,'defaulttextinterpreter','none');
            step_len = 360/parms.num_of_direction_bins;
            y_presentation_axis=-180:step_len:180;
            bar(y_presentation_axis,Cell.MD.rate_hist);
            xlabel('angle');
            ylabel('count');
            strTitle = sprintf('Moving_Directionality ray=%0.3f',Cell.MD.rayleigh.score);
            title(strTitle);
            
            subplot(n,m,2);
            imagesc(rate_mat);
            axis image
            hold on;
            
            max_inds_len =length(max_inds);
            
            for cen = 1:max_inds_len;
                plot((max_inds(cen,2)),(max_inds(cen,1)), 'x', 'color','k','MarkerSize',12)
            end
            
            
            
            subplot(n,m,3)
            imagesc(zone_mat);
            axis image;
            colorbar;
            strTitle = sprintf('rate diff b/w borders and nonborder= %0.3f \n top =%0.3f, bottom=%0.3f, left=%0.3f, right=%0.3f', S.border_rate_diff, S.top, S.bottom, S.left, S.right);
            title(strTitle);
            
            subplot(n,m,4);
            plot(Cell.sorted_means, 'o-')
            ymax= max(S.sorted_means);
            xmax= length(S.sorted_means)+1;
            axis ([0 xmax 0 ymax])
            hold on;
            %            y = S.bkgd_firing;
            %            line('XData', [0 xmax], 'YData', [y y], 'LineStyle', '-.', ...
            %            'LineWidth', 2, 'Color','r')
            strTitle= sprintf('avg mean=%0.3f \n variability = %0.3f \n mean norm mean=%0.3f \n fano factor= %0.3f', Cell.overall_mean, Cell.variance, mean(Cell.norm_means), S.var_over_mean);
            title(strTitle);
            
            six_orientation_pts = find_six_points(rate_mat_after_auto_pearson, Cell.PF_radius)
            Cell.six_orientation_pts = six_orientation_pts;
            
            fig_pair= subplot(n,m,5);
            fig_pair= plot_ellipse(rate_mat_after_auto_pearson, Cell.module.major, Cell.module.minor ,Cell.module.phi, six_orientation_pts, fig_pair);
            
            
            
            draw_orientation_line= subplot(n,m,6);
            
            %%%COMMENT BACK IN AFTERWARDS.......................
            
            [draw_orientation_line,smallest_degree, smallest_degree_wall] = orientation_analysis(rate_mat_after_auto_pearson, rate_mat, Cell.six_orientation_pts, draw_orientation_line, max_inds) ;
            
            
            
            Cell.smallest_degree= smallest_degree;
            Cell.smallest_degree_wall= smallest_degree_wall;
            
            
            
            %no_angle_change_now=0;
            arena_phi=0;
            
            [min_angle, main_axis, axis_and_orient] = adjusted_orientation_analysis(rate_mat_after_auto_pearson, rate_mat, Cell.six_orientation_pts, arena_phi ) ;
            
            Cell.min_angle=min_angle; %adjusted angle
            Cell.main_axis=main_axis; %should be same as adjusted angle, calculated differently
            Cell.axis_and_orient=axis_and_orient;
            
            location_inds = find(S.location == 1);
            max_peak_location=location_inds;
            
            subplot(n,m,6)
            strTitle= sprintf('smallest degree = %0.3f \n wall =%d \n location=%d %d=', ...
                smallest_degree, smallest_degree_wall, max_peak_location);
            title(strTitle);
            
            
            %find head direction
            
            %%remove after- only for no arena change
            
            pos_y2= Cell.pos.y2;
            pos_x2=Cell.pos.x2;
            pos_y= Cell.pos.y;
            pos_x=Cell.pos.x;
            
            pos_HD =wrapToPi(atan2(pos_y2-pos_y,pos_x2-pos_x));
            
            %% correct HD
            [location_of_gaze_x, location_of_gaze_y, hist_count_gaze, gaze_rate, gaze_spike_rate]=findLocationGaze(pos_HD,parms,rate_mat,Cell.pos.t, pos_mean_y, pos_mean_x, S.spk.t);
            
            Cell.location_of_gaze_x=location_of_gaze_x;
            Cell.location_of_gaze_y=location_of_gaze_y;
            Cell.hist_count_gaze= hist_count_gaze;
            
            gaze_spike_rate= gaze_spike_rate./gaze_rate;  %spike rate dividied by amount of times gazed
            
            gaze_rate= SmoothRateMat(gaze_rate,parms);
            
            gaze_spike_rate= SmoothRateMat(gaze_spike_rate,parms);
            
            Cell.gaze_spike_rate=gaze_spike_rate;
            Cell.gaze_rate= gaze_rate;
            
            
            subplot(n,m,7)
            imagesc(gaze_rate);
            axis equal; axis off;
            title('gaze rate map');
            
            subplot(n,m,8)
            imagesc(gaze_spike_rate);
            axis equal; axis off;
            title('gaze rate spike map');
            
            
            subplot(n,m,9)
            imagesc(zone_spike_mat);
            axis equal; axis off;
            title('zone mat by number (not rate) of spikes');
            
            
            subplot(n,m,9)
            imagesc(zone_mat);
            axis equal; axis off;
            
            % divide peak firing rate by area of place field
            count=1;
            
            for h=1:length(S.sorted_means)
                PF_area= sum(S.zone_mat==S.sorted_means(h));
                PF_area=sum(PF_area);
                PF_areas(count)= PF_area; % areas of place fields in same order as sorted means
                count=count+1;
            end
            
            max_PF_area= max(PF_areas); %maximum PF size to be used at standard PF size
            rate_over_area= S.sorted_means.*(PF_areas/max_PF_area); %multiple rates by percentage of PF size
            
            subplot(n,m,10)
            plot(rate_over_area, 'o-')
            ymax= max(rate_over_area);
            xmax= length(rate_over_area)+1;
            axis ([0 xmax 0 ymax])
            hold on;
            title('sorted firing rates multiplied by percentage of PF size')
            
            Cell.rate_over_area=rate_over_area;
            Cell.PF_areas= PF_areas;
            
            rate_over_area_mat = S.zone_mat;
            
            for h=1:length(rate_over_area);
                rate_over_area_mat(rate_over_area_mat==S.sorted_means(h)) = rate_over_area(h);
            end
            
            Cell.rate_over_area_mat=rate_over_area_mat;
            
            [PF_sorted_by_size, SortIndex] = sort(PF_areas);
            firing_rates_by_PF_size = S.sorted_means(SortIndex);
            rate_over_area_by_PF_size= rate_over_area(SortIndex);
            
            Cell.firing_rates_by_PF_size= firing_rates_by_PF_size;
            Cell.rate_over_area_by_PF_size=rate_over_area_by_PF_size;
            
            subplot(n,m,7)
            imagesc(rate_over_area_mat);
            axis equal; axis off;
            
            subplot(n,m,8)
            scatter(PF_areas, S.sorted_means)
            lsline;
            title('firing rates vs. PF size')
            
            subplot(n,m,9)
            scatter(PF_areas, rate_over_area)
            lsline;
            title('firing rates over area vs. PF size')
            
            disp('');
            
            
            
            cd(parms.dir_save_pictures);
            %  saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %
            
            close all;
            
        end
        
        S=Cell;
        
        bSave = false;
        if (bSave)
            cd(parms.dir_save_data);
            save(sprintf('Cell_r%s_d%s_s%s_t%d_c%d.mat',Cell.rat,Cell.date,...
                Cell.session,Cell.tetrode,Cell.cell),'S');
            %     %         % save pictuers % debugger - return
            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
        end
        
        
        count=count+1;
        % alpha
        %      close all;
    end
    
    % saves results
    cd(parms.dir_save_data);
    title1 =strcat (file_name);
    %title1 =strcat ('Results_',file_name);
    save (title1,'S');
    
end


count = count - 1;


    function S= get_zones(rate_mat, S)
        %%find max points in rate mat
        %%find max points in auto mat, use this to find radius of PF
        %%use PF radius in rate mat to create zone mat by comparing max pts to PF
        %%radius
        
        % calculate maximum points in rate mat
        
        rate_mat_orig = rate_mat;
        rate_mat(isnan(rate_mat))=0;    %change nans to zeros
        
        [size_x,size_y]=size(rate_mat);
        
        % embed rate_mat into a matrix with zeros on edges
        %[so points on edges can be found]
        
        new_rate_mat = zeros(size_x+2,size_y+2);
        
        new_rate_mat(2:size_x+1,2:size_y+1)=rate_mat;
        
        rate_mat = new_rate_mat;
        [size_x,size_y]=size(rate_mat);
        
        % look for local maxima in rate_mat
        
        max_inds_len = 0;
        max_inds = [];
        
        for fig_i = 2:size_x-1
            for j = 2:size_y-1
                if rate_mat(fig_i,j) > rate_mat(fig_i+1,j) && ...
                        rate_mat(fig_i,j) > rate_mat(fig_i-1,j) && ...
                        rate_mat(fig_i,j) > rate_mat(fig_i,j+1) && ...
                        rate_mat(fig_i,j) > rate_mat(fig_i,j-1)
                    hold on
                    max_inds_len = max_inds_len+1;
                    max_inds(max_inds_len,1) = fig_i-1;     %list of maximum pt indices
                    max_inds(max_inds_len,2) = j-1;
                end
            end
        end
        
        rate_mat = rate_mat_orig;       %return orignal rate mat without zeros border
        [size_x,size_y]=size(rate_mat);
        
        
        
        
        
        %create auto correlated rate map
        
        rate_mat_after_auto_pearson = Cross_Correlation(rate_mat,rate_mat);
        Cell.rate_mat_after_auto_pearson = rate_mat_after_auto_pearson;
        
        
        %find maximum points in rate_mat_after_auto_pearson
        
        %change autocorr to redone autocorrelation
        Cell.autocorr = rate_mat_after_auto_pearson;
        
        [size_x, size_y] = size(Cell.autocorr);
        auto_max_inds_len = 0;
        
        %find autocorrelation mat maximum points
        
        for fig_i = 2:size_x-1
            for j = 2:size_y-1
                if Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i+1,j) && ...
                        Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i-1,j) && ...
                        Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i,j+1) && ...
                        Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i,j-1)
                    % plot(j,fig_i,'x');
                    
                    auto_max_inds_len = auto_max_inds_len+1;
                    auto_max_inds(auto_max_inds_len,1) = fig_i;   %indices of maximum pts
                    auto_max_inds(auto_max_inds_len,2) = j;
                end
            end
        end
        
        S.rate_mat_after_auto_pearson = rate_mat_after_auto_pearson;
        
        S.auto_max_inds= auto_max_inds;
        
        % find radius of PF using half distance between center and closest max pt
        
        cen = 1:auto_max_inds_len;
        distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
        
        new_distances = sort(distances (:));
        PF_radius = new_distances(2)/2; %% finds 2nd min since min is 0
        
        % has to be at least 70% of calculated radius distance
        
        PF_radius = PF_radius; %% should be correct.
        
        S.PF_radius = PF_radius;
        
        
        % find distances that are too close together between peaks
        
        h=1;
        too_close=[];
        
        for cen = 1:max_inds_len-1
            for cen2= (cen+1):max_inds_len
                peak_distance = Distance(max_inds(cen,1), max_inds(cen,2),max_inds(cen2,1), max_inds(cen2,2));
                if peak_distance < PF_radius
                    too_close (h,1) = cen;
                    too_close(h,2) = cen2;
                    
                    h=h+1;
                end
            end
        end
        
        if ~isempty (too_close)
            
            [size_x, size_y] = size(too_close);
            
            too_close_len = size_x;
            
            
            %remove one of the peaks that are too close with lower firing rate
            
            h=1;
            
            for cen = 1: too_close_len
                if rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) >=...
                        rate_mat((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ;
                    
                    remove(h) = too_close(cen,2);
                    
                    h= h+1;
                    
                elseif rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) <...
                        rate_mat((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ;
                    
                    remove(h)= too_close(cen,1);
                    
                    h=h+1;
                    
                else
                    disp ('wtf')
                end
            end
            
            max_inds(remove,1) = 0;
            max_inds(remove,2) = 0;
            
            max_inds(all(max_inds==0,2),:) = [];
            
        end
        
        
        % give each max point a value
        
        [size_x, size_y] = size(rate_mat);
        
        max_inds_len= length(max_inds);
        
        zone_mat = zeros(size_x, size_y);
        
        zone_num = 1;
        
        for cen = 1:max_inds_len;
            zone_mat((max_inds(cen,1)),(max_inds(cen,2))) = zone_num;
            
            zone_num = zone_num+1;
        end
        
        
        peak_zone_mat = zone_mat;
        
        % if distance to max pt is less than 70% PF radius, assign it value at max pt
        
        PF_radius = 0.8* PF_radius;
        
        [size_x, size_y] = size(rate_mat);
        
        for cen=1:max_inds_len;
            for fig_i =1:size_x
                for j =1:size_y
                    if Distance(fig_i, j, max_inds(cen,1), max_inds(cen,2)) < PF_radius
                        zone_mat(fig_i,j)= cen;
                    end
                end
            end
        end
        
        number_zone_mat = zone_mat; %zone_mat with arbituary numbers instead of firing rates
        
        number_zones = (unique(peak_zone_mat))' ;  %number of PF and bkgd, includes 0
        
        number_zones(number_zones==0)= [];
        
        number_zones_len = length(number_zones);
        
        %
        % %find mean firing rate for each zone
        %% finds mean firing rates of each zone
        
        % rate_mat(isnan(rate_mat))= 0;
        %
        % mean_values_len = 0;
        %
        % for cen = 1:number_zones_len;
        %
        %       find_spk = find(number_zone_mat == number_zones(cen));
        %       if sum(find_spk>length(rate_mat))>0
        %         disp('')
        %       end
        %
        %       mean_values_len= mean_values_len+1;
        %       mean_values_list(mean_values_len) = mean(rate_mat(find_spk));
        
        %% finds peak firing rate of each grid cell place field
        
        rate_mat(isnan(rate_mat)) = 0;
        
        peak_rates_len = 0;
        
        peak_rates_list = [];
        
        for cen = 1:max_inds_len
            find_spk = find(peak_zone_mat == cen);
            if sum(find_spk>length(rate_mat))>0
                disp('')
            end
            
            peak_rates_len = peak_rates_len+1;
            peak_rates_list(peak_rates_len)= rate_mat(find_spk);
            
            % orig_mean_values_list = mean_values_list; %use orig for color coding plot by firing rate
            
            % mean_values_list(mean_values_list==0) = []; %remove 0s from removed too-small PFs
            % mean_values_list(isnan(mean_values_list)) = []; %remove NaNs
        end
        
        
        %%changes zone map to peak firing rate values
        
        for cen=1:number_zones_len;
            zone_mat(find(zone_mat==number_zones(cen))) = peak_rates_list(cen);
        end
        
        
        %         % zone mat based on number of spikes (not firing rate)
        %         % uses spike_mat instead of rate_mat
        %
        %         zone_spike_mat= number_zone_mat;
        %
        %          for cen=1:number_zones_len;
        %             spike_sum= sum(spike_mat(find(number_zone_mat==cen)));
        %             zone_spike_mat(find(number_zone_mat==number_zones(cen))) = spike_sum;
        %          end
        %
        %         S.zone_spike_mat=zone_spike_mat;
        
        S.max_inds=max_inds;
        
        %bkgd_firing = mean_values_list(1);
        
        %S.bkgd_firing = bkgd_firing;
        
        %mean_values_list(1) = []; %drop background mean firing rate
        
        S.number_of_PF = length(peak_rates_list);
        
        S.overall_mean = mean(peak_rates_list); %overall average mean firing rate of cell
        
        S.mean_values_list= peak_rates_list';
        %
        S.variance = var(peak_rates_list);
        
        PF_sum= nnz(number_zone_mat);
        
        PF_coverage = PF_sum/(size_x*size_y);
        
        S.PF_coverage = PF_coverage;
        
        S.zone_mat = zone_mat;
        
        % plot mean values
        sorted_means = sort(peak_rates_list);
        
        S.sorted_means=sorted_means;
        
        norm_means = S.sorted_means/max(S.sorted_means);
        
        S.norm_means = norm_means;
        
        S.variance_norm_means = var(norm_means);
        
        %storing vars back to S:
        S.new_rate_mat=rate_mat;
        S.size_x=size_x;
        S.size_y=size_y;
        
        var_over_mean = var(sorted_means)/mean(sorted_means);
        S.var_over_mean = var_over_mean;
        
        %%
        
        top = unique(number_zone_mat(1,:));
        top(find(top==0))= [];
        
        left = unique(number_zone_mat(:,1));
        left(find(left==0))= [];
        
        right = unique(number_zone_mat(:,size_y));
        right(find(right==0))= [];
        
        bottom = unique(number_zone_mat(size_x,:)) ;
        bottom(find(bottom==0))= [];
        
        A = union(top, left);
        C = union(bottom, right);
        
        border_zones = union(A, C);
        border_zones(find(border_zones==0))= [];
        
        border_mean_rates = peak_rates_list(border_zones);
        
        nonborder_means = setdiff(peak_rates_list, border_mean_rates);
        
        border_to_nonborder_diff= mean(border_mean_rates)/mean(nonborder_means);
        
        S.border_rate_diff = border_to_nonborder_diff;
        
        %% looks at border to nonborder difference with max peak removed
        
        peak_rate = max(peak_rates_list);
        
        if any(border_mean_rates==peak_rate)
            nonborder_means_wo_max=nonborder_means;
            border_means_wo_max= border_mean_rates;
            border_means_wo_max(border_means_wo_max==peak_rate)=[];
        elseif any(nonborder_means ==peak_rate)
            nonborder_means_wo_max=nonborder_means;
            border_means_wo_max= border_mean_rates;
            nonborder_means_wo_max(nonborder_means_wo_max==peak_rate)=[];
        else
            disp('ERROR ERROR')
        end
        
        border_to_nonborder_diff_wo_max= mean(border_means_wo_max)/mean(nonborder_means_wo_max);
        
        S.border_to_nonborder_diff_wo_max= border_to_nonborder_diff_wo_max;
        
        
        top_to_central = mean(peak_rates_list(top))/mean(nonborder_means);
        S.top = top_to_central;
        
        bottom_to_central = mean(peak_rates_list(bottom))/mean(nonborder_means);
        S.bottom = bottom_to_central;
        
        left_to_central= mean(peak_rates_list(left))/mean(nonborder_means);
        S.left= left_to_central;
        
        right_to_central= mean(peak_rates_list(right))/mean(nonborder_means);
        S.right = right_to_central;
        
        %% finds position of highest firing rate zone
        
        zone_num = find(peak_rates_list==max(peak_rates_list));
        
        location(1) = any(left==zone_num);
        location(2) = any(right==zone_num);
        location(3) = any(bottom==zone_num);
        location(4) = any(top==zone_num);
        
        where = sum(location);
        
        S.location = location;
        
        S.where = where;
        
        S.peak_rates_list = peak_rates_list;
        
        
        max_firing_index(1,1) = max_inds(zone_num, 1);
        max_firing_index(1,2) = max_inds(zone_num, 2);
        
        S.max_index = max_firing_index;
        
        [size_x, size_y] = size(rate_mat);
        
        norm_max_index(1,1) = max_firing_index(1,1) / size_x;
        
        norm_max_index(1,2) = max_firing_index(1,2) / size_y;
        
        S.norm_max_index = norm_max_index;
        
        S.zone_mat = zone_mat;
        
        S.number_zone_mat = number_zone_mat ;
        
        S.peak_zone_mat = peak_zone_mat;
        
        S.i = i;
        
        %% finds distance of max peak from border
        
        h = 1;
        
        for cen = [1, size_x]
            for cen2 = 1:size_y
                max_peak_distance_top_bottom (1,h) = Distance(max_firing_index(1, 1),max_firing_index(1, 2), cen , cen2);
                
                h = h+1;
            end
        end
        
        h = 1;
        
        for cen = [1, size_y]
            for cen2= 1:size_x
                max_peak_distance_right_left(1,h)= Distance(max_firing_index(1, 1),max_firing_index(1, 2), cen , cen2);
                
                h= h+1;
            end
        end
        
        max_peak_distance = union(max_peak_distance_top_bottom, max_peak_distance_right_left);
        
        max_peak_distance = min(max_peak_distance);
        
        S.max_peak_distance = max_peak_distance;
        
        disp('');
        %% finds the points for the module ellipse
        
        [size_x, size_y] = size(S.autocorr);
        
        auto_distances = [];
        new_distances = [];
        
        %%finds distances from center of auto_corr_map to all max peaks
        
        cen = 1:length(S.auto_max_inds);
        auto_distances = Distance(S.auto_max_inds(cen, 1),S.auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
        
        new_distances = sort(auto_distances (:));
        
        
        disp('');
        
    end %of getzones

disp('');

end % of main script












