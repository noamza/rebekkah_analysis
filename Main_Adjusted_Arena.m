function Main
% Date: 24 of July 2014
dbstop if error

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\results gridness 0.3 with HD';
parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\all examples 22_03_15 images';
% parms.dir_save_pictures2= 'C:\Users\Dori\Desktop\Rebekkah data\results testing';
% parms.dir_save_pictures3= 'C:\Users\Dori\Desktop\Rebekkah data\results testing';
parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\all examples 22_03_15';

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
    
    if ~isempty(Cell.pos.x2)
        dt = Cell.pos.t(2)-Cell.pos.t(1);
        
        % calculate the the rat's head direction (using average of both leds)
        pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
        % get rate matrix - using function: Creat_Rate_Map
        rate_mat=CreateRateMap(Cell.pos.x,Cell.pos.y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        Cell.rate_mat = rate_mat;
        
        S=get_zones(rate_mat,dat.S);
        
        % use get matrix after auto pearson correlation - using function: Auto_Pearson_Correlation
        rate_mat_after_auto_pearson = AutoPearsonCorrelation(rate_mat);
        
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
        Cell.six_orientation_pts= S.six_orientation_pts;
        Cell.max_inds_adjusted=S.max_inds_adjusted;
        
        
        Cell_num(count)=i;
        
        % plots and graphs figures
        
        
        bPrintPicture = true;
        if (bPrintPicture)
            % plot all the picturs and graphs
            strTitle = sprintf('Cell_i=%d_r%s_d%s_s%s_t%d_c%d',i,Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell);
            fig =figure('name',strTitle);
            % title(strTitle);
            % plot the rat's path (and the spikes on it)
            n = 3;
            m = 3;
            subplot(n,m,1);
            plot(pos_mean_x,pos_mean_y,'k');hold on;
            plot(spk_x,spk_y,'.r');
            axis equal;axis off;
            axis ij
            title(file_name);
            
            hold on;
            
            
            
            %             subplot(n,m,2);
            %             imagesc(rate_mat); axis image
            %
            %             subplot(n,m,3);
            %             imagesc(rate_mat_after_auto_pearson); axis image
            %             hold on;
            %             title(sprintf('gridness=%0.2f /n elis_grid=%0.2f', S.gridness2, gridness3));
            
            %             subplot(n,m,4);
            %             hist(Hist_Moving_Directionality_In_Head_Directionality,ax);
            %             title('Moving Directionality In Head Directionality');
            
            %             subplot(n,m,4);
            %             set(0,'defaulttextinterpreter','none');
            %             step_len = 360/parms.num_of_direction_bins;
            %             y_presentation_axis=-180:step_len:180;
            %             bar(y_presentation_axis,Cell.HD.rate_hist);
            %             xlabel('angle');
            %             ylabel('count');
            %             strTitle = sprintf('Head_Directionality ray=%0.3f',Cell.HD.rayleigh.score);
            %             title(strTitle);
            
            %             subplot(n,m,6);
            %             set(0,'defaulttextinterpreter','none');
            %             step_len = 360/parms.num_of_direction_bins;
            %             y_presentation_axis=-180:step_len:180;
            %             bar(y_presentation_axis,Cell.MD.rate_hist);
            %             xlabel('angle');
            %             ylabel('count');
            %             strTitle = sprintf('Moving_Directionality ray=%0.3f',Cell.MD.rayleigh.score);
            %             title(strTitle);
            
            subplot(n,m,2);
            imagesc(rate_mat);
            axis image
            hold on;
            
            max_inds_len =length(S.max_inds);
            
            for cen = 1:max_inds_len;
                plot((max_inds(cen,2)),(max_inds(cen,1)), 'x', 'color','k','MarkerSize',12)
            end
            
            
            
            subplot(n,m,3)
            imagesc(zone_mat);
            axis image;
            colorbar;
            strTitle = sprintf('# of PF=%0.3d \n radius len =%0.3f \n PF coverage=%0.3f \n rate difference b/w borders and nonborder= %0.3f \n top =%0.3f, bottom=%0.3f, left=%0.3f, right=%0.3f', Cell.number_of_PF, Cell.PF_radius, Cell.PF_coverage, S.border_rate_diff, S.top, S.bottom, S.left, S.right);
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
            
            six_orientation_pts = find_six_points(Cell.autocorr, Cell.PF_radius)
            Cell.six_orientation_pts = six_orientation_pts;
            
            fig_pair= subplot(n,m,5);
            fig_pair= plot_ellipse(Cell.autocorr, Cell.module.major, Cell.module.minor ,Cell.module.phi, Cell.six_orientation_pts, fig_pair);
            
            
            
            draw_orientation_line= subplot(n,m,6);
            [draw_orientation_line,smallest_degree, smallest_degree_wall] = orientation_analysis(Cell.autocorr, rate_mat, Cell.six_orientation_pts, draw_orientation_line, S.max_inds) ;
            
            [min_top, min_bottom] = orientation_analysis(Cell.autocorr, rate_mat, Cell.six_orientation_pts, draw_orientation_line, S.max_inds) ;
            
            
            Cell.smallest_degree= smallest_degree;
            Cell.smallest_degree_wall= smallest_degree_wall;
            Cell.min_top= min_top;
            Cell.min_bottom= min_bottom;
            
            
            
            %% plot fitted arena around trajectory map
            
            subplot(n,m,1)
            [rotated_sq_x, rotated_sq_y, new_angle, adjusted_rate_mat] = draw_fitted_arena_box(file_name, parms.dir_load_data, Cell.smallest_degree, Cell.rate_mat);
            
            x_1_r=rotated_sq_x(1);
            x_2_r=rotated_sq_x(2);
            x_3_r=rotated_sq_x(3);
            x_4_r=rotated_sq_x(4);
            
            y_1_r=rotated_sq_y(1);
            y_2_r=rotated_sq_y(2);
            y_3_r=rotated_sq_y(3);
            y_4_r=rotated_sq_y(4);
            
            line([x_1_r, x_3_r],[y_1_r, y_3_r])
            line([x_3_r, x_4_r],[y_3_r, y_4_r])
            line([x_4_r, x_2_r], [y_4_r, y_2_r])
            line([x_2_r, x_1_r], [y_2_r, y_1_r])
            
            subplot(n,m,2)
            
            line([x_1_r, x_3_r],[y_1_r, y_3_r])
            line([x_3_r, x_4_r],[y_3_r, y_4_r])
            line([x_4_r, x_2_r], [y_4_r, y_2_r])
            line([x_2_r, x_1_r], [y_2_r, y_1_r])
            
            Cell.adjusted_rate_mat = adjusted_rate_mat;
           
            autocorr_adjusted= AutoPearsonCorrelation(adjusted_rate_mat);
            
            draw_orientation_line= subplot(n,m,7);
            [draw_orientation_line,smallest_degree_adjusted, smallest_degree_adjusted, smallest_degree_wall_adj] = ...
                 orientation_analysis(autocorr_adjusted_angle, adjusted_rate_mat, Cell.six_orientation_pts, draw_orientation_line, S.max_inds_adjusted) ;
            
             Cell.smallest_degree_wall_adjusted= smallest_degree_wall_adjusted;
             
             subplot(n,m,8)
             imagesc(autocorr_adjusted_angle);
            
            
            disp('');
            
            %%%%% RECOMMENT AFTER.
            
            %            if var(S.sorted_means)/mean(S.sorted_means) < 1
            %            cd(parms.dir_save_pictures);
            %            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
            %            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %
            %
            %            elseif var(S.sorted_means)/mean(S.sorted_means) > 1 &&  var(S.sorted_means)/mean(S.sorted_means) < 2
            %            cd(parms.dir_save_pictures2);
            %            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
            %            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %
            %
            %            elseif var(S.sorted_means)/mean(S.sorted_means) > 2
            %            cd(parms.dir_save_pictures3);
            %             saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
            %            saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.jpg',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %
            %
            %             else
            %                 disp('wtf');
            %             end
            
            
            cd(parms.dir_save_pictures);
            %   saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
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




    function S= get_zones(rate_mat, adjusted_rate_mat, S)
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
        
        
        [size_x,size_y]=size(adjusted_rate_mat);
        
        % look for local maxima in adjusted rate_mat
        
        max_inds_len_adj = 0;
        max_inds_adjusted = [];
        
        for fig_i = 2:size_x-1
            for j = 2:size_y-1
                if adjusted_rate_mat(fig_i,j) > adjusted_rate_mat(fig_i+1,j) && ...
                        adjusted_rate_mat(fig_i,j) > adjusted_rate_mat(fig_i-1,j) && ...
                        adjusted_rate_mat(fig_i,j) > adjusted_rate_mat(fig_i,j+1) && ...
                        adjusted_rate_mat(fig_i,j) > adjusted_rate_mat(fig_i,j-1)
                    hold on
                    max_inds_len_adj = max_inds_len_adj+1;
                    max_inds_adjusted(max_inds_len_adj,1) = fig_i-1;     %list of maximum pt indices
                    max_inds_adjusted(max_inds_len_adj,2) = j-1;
                end
            end
        end
        
        S.max_inds_adjusted=max_inds_adjusted;
        
        %create auto correlated rate map
        
        rate_mat_after_auto_pearson = AutoPearsonCorrelation(adjusted_rate_mat);
        Cell.rate_mat_after_auto_pearson = rate_mat_after_auto_pearson;
        
        
        %find maximum points in rate_mat_after_auto_pearson
        
        [size_x, size_y] = size(Cell.autocorr);
        auto_max_inds_len = 0;
        
        %find autocorrelation mat maximum points
        
        for fig_i = 2:size_x-1
            for j = 2:size_y-1
                if rate_mat_after_auto_pearson(fig_i,j) > rate_mat_after_auto_pearson(fig_i+1,j) && ...
                        rate_mat_after_auto_pearson(fig_i,j) > rate_mat_after_auto_pearson(fig_i-1,j) && ...
                        rate_mat_after_auto_pearson(fig_i,j) > rate_mat_after_auto_pearson(fig_i,j+1) && ...
                        rate_mat_after_auto_pearson(fig_i,j) > rate_mat_after_auto_pearson(fig_i,j-1)
                    % plot(j,fig_i,'x');
                    
                    auto_max_inds_len = auto_max_inds_len+1;
                    auto_max_inds_adj(auto_max_inds_len,1) = fig_i;   %indices of maximum pts
                    auto_max_inds_adj(auto_max_inds_len,2) = j;
                end
            end
        end
        
        
        
        S.rate_mat_after_auto_pearson = rate_mat_after_auto_pearson;
        
        S.auto_max_inds_adj= auto_max_inds_adj;
        
        
        
        % find radius of PF using half distance between center and closest max pt
        
        cen = 1:auto_max_inds_len;
        distances = Distance(auto_max_inds_adj(cen, 1),auto_max_inds_adj(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
        
        new_distances = sort(distances (:));
        PF_radius = new_distances(2)/2; %% finds 2nd min since min is 0
        
        % has to be at least 70% of calculated radius distance
        
        PF_radius = PF_radius;
        
        S.PF_radius = PF_radius;
        
        
        % find distances that are too close together between peaks
        
        h=1;
        too_close=[];
        
        for cen = 1:max_inds_len-1
            for cen2= (cen+1):max_inds_len
                peak_distance = Distance(max_inds(cen,1), max_inds(cen,2),max_inds(cen2,1), max_inds(cen2,2));
                if peak_distance < 1.25*(PF_radius)
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
                if adjusted_rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) >...
                        adjusted_rate_mat((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ;
                    
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
        
        PF_radius = 0.65* PF_radius;
        
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
        
        location1 = any(top==zone_num);
        location2 = any(bottom==zone_num);
        location3 = any(left==zone_num);
        location4 = any(right==zone_num);
        
        location(1)= location1;
        location(2)= location2;
        location(3) = location3;
        location(4) = location4;
        
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
        %%% finds the points for the module ellipse
        
        [size_x, size_y] = size(S.autocorr);
        
        auto_distances = [];
        new_distances = [];
        
        %%finds distances from center of auto_corr_map to all max peaks
        
        cen = 1:length(S.auto_max_inds);
        auto_distances = Distance(S.auto_max_inds(cen, 1),S.auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
        
        new_distances = sort(auto_distances (:));
        
        % finds the 6 closest pts to the center
        
        auto_dist_inds1 = [];
        auto_dist_inds2 = [];
        auto_dist_inds3 = [];
        auto_dist_inds4 = [];
        auto_dist_inds5 = [];
        auto_dist_inds6 = [];
        
        auto_dist_inds1 = find(auto_distances == new_distances(2));
        auto_dist_inds2 = find(auto_distances == new_distances(3));
        auto_dist_inds3 = find(auto_distances == new_distances(4));
        auto_dist_inds4 = find(auto_distances == new_distances(5));
        auto_dist_inds5 = find(auto_distances == new_distances(6));
        auto_dist_inds6 = find(auto_distances == new_distances(7));
        
        union1 = union(auto_dist_inds1, auto_dist_inds2);
        union2 = union(auto_dist_inds3, auto_dist_inds4);
        union3 = union(auto_dist_inds5, auto_dist_inds6);
        union1 = union(union1, union2);
        
        auto_dist_inds = union(union1, union3);
        
        for k= 1:length(auto_dist_inds);
            six_orientation_pts (k,1) = S.auto_max_inds(auto_dist_inds(k),1);
            six_orientation_pts (k,2) = S.auto_max_inds(auto_dist_inds(k),2);
        end
        
        %adds center point
        
        six_orientation_pts (7,1) = (size_x/2)+0.5;
        six_orientation_pts (7,2) = (size_y/2)+0.5;
        
        S.six_orientation_pts= six_orientation_pts;
        
        
        
        
    end %of getzones



disp('');

end % of main script












