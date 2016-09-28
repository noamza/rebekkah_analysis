function Main
% Date: 24 of July 2014
dbstop if error
 
parms.dir_load_data = 'C:\Users\Dori\Desktop\elipse plots PF 8 results'; 
parms.dir_save_pictures='N:\users\rebekkah\final data smoothed\examine fano factor first half second half\images';
parms.dir_save_data = 'N:\users\rebekkah\final data smoothed\examine fano factor first half second half\results'; 
  
dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
 
parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120; %midpoint1 (+pi/2),midpoint2(-pi/2)
parms.bin_size=3;
parms.sigma = 1.5;

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
    
    begin_inds=[];
    end_inds=[];
    pos_t_len=[];
    pos_t_end=[];
    pos_t_begin=[];
    
    pos_t_len= length(Cell.pos.t);
    
    if mod(pos_t_len,2) ==0   % Even  number
        
        begin_inds= 1:1:pos_t_len/2;
        end_inds= (pos_t_len/2)+1:1:pos_t_len;
    elseif mod(pos_t_len,2) == 1    % odd number
        
        begin_inds= 1:1:(pos_t_len/2)+0.5;
        end_inds= (pos_t_len/2)+1.5:1:pos_t_len;
    end
    
    
    pos_t_begin= Cell.pos.t(begin_inds);   %first half of session
    pos_t_end= Cell.pos.t(end_inds); %second half of session
   
    tmp_begin = abs(Cell.spk.t-pos_t_begin(end)); %subtract last timestamp of first session
    tmp_end= abs(Cell.spk.t- pos_t_end(1)) ;% subtract first timestamp of second session 
    [~, idx_begin] = min(tmp_begin); %index of first session
    [~, idx_end] = min(tmp_end) ;%index of first session
    spk_t_begin = Cell.spk.t(1:idx_begin); %spike vector of first sesh
    spk_t_end = Cell.spk.t(idx_end+1:length(Cell.spk.t)); %spike vector of second sesh
    
    if ~isempty(Cell.pos.x2)
        dt = Cell.pos.t(2)-Cell.pos.t(1);
        
        % calculate the the rat's head direction (using average of both leds)
        pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t); 
        
        % divide into first 5 minutes and last 5 minutes 
        pos_mean_x_begin= pos_mean_x(begin_inds); %x position for first half of session
        pos_mean_x_end= pos_mean_x(end_inds); %x position for second half of session
        pos_mean_y_begin= pos_mean_y(begin_inds); %y position for first half of session
        pos_mean_y_end= pos_mean_y(end_inds); %y position for second half of session
        
        % confirm the Cell.spk.t and spk_x are same size and correlate to
        % eachother!!!!!!!!!!!!!!!
        
        spk_x_begin= spk_x(1:idx_begin);
        spk_x_end= spk_x(idx_end:length(spk_x));
        spk_y_begin= spk_y(1:idx_begin);
        spk_y_end= spk_y(idx_end:length(spk_y));
       
        % get rate matrix - using function: Creat_Rate_Map 
        rate_mat_begin=CreateRateMap(pos_mean_x_begin,pos_mean_y_begin,pos_t_begin,spk_x_begin,spk_y_begin,spk_t_begin,parms);
        rate_mat_end=CreateRateMap(pos_mean_x_end,pos_mean_y_end,pos_t_end,spk_x_end,spk_y_end,spk_t_end,parms);
        rate_mat_total=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        Cell.rate_mat_begin = rate_mat_begin;
        Cell.rate_mat_end = rate_mat_end;
        Cell.rate_mat_total= rate_mat_total;
        
        S_begin=get_zones(rate_mat_begin, dat.S);
        S_end=get_zones(rate_mat_end, dat.S);
        S_total=get_zones(rate_mat_total, dat.S);
       
%         % use get matrix after auto pearson correlation - using function: Auto_Pearson_Correlation
%         
%         rate_mat_after_auto_pearson_begin = AutoPearsonCorrelation(rate_mat_begin);
%         rate_mat_after_auto_pearson_end = AutoPearsonCorrelation(rate_mat_end);
%         
%         Cell.rate_mat_after_auto_pearson_begin = rate_mat_after_auto_pearson_begin;
%         Cell.rate_mat_after_auto_pearson_end = rate_mat_after_auto_pearson_end;
%         
        % don't use this gridness score, so removed.
        
        % R_outer = FindROuter(rate_mat_after_auto_pearson,parms);
        % TODO: to print value of gridness2
        % TODO: to check if the figure print is really required
        % [gridness3] = GridnessRadius(rate_mat_after_auto_pearson,parms,R_outer,i);
        
        % Cell.gridness2 = gridness3;
        % pop_gridness2(count)=gridness3;
        
        % calculate hist of HD
       % Hist_Moving_Directionality_In_Head_Directionality = ComputePrecentOfMovingDirectionalityInHeadDirectionality(Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2);
        %ax = deg2rad(-180:180);
        %hist(Hist_Moving_Directionality_In_Head_Directionality,ax);
       % [~, location] = max(hist(Hist_Moving_Directionality_In_Head_Directionality,...
%             ax));
%         alpha = ax(location);
%         Cell.alpha = alpha;
        
        % calculate HD
%         [HD_rate_hist,HD_rayleigh_score,HD_rayleigh_angle]=ComputeHeadDirectionality...
%             (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2,Cell.spk.t,parms);
%         pop_rayleighHD_score(count)=HD_rayleigh_score;
%         pop_rayleighHD_angle(count)=HD_rayleigh_angle;
%         Cell.HD.rate_hist = HD_rate_hist;
%         Cell.HD.rayleigh.score=HD_rayleigh_score;
%         Cell.HD.rayleigh.angle=HD_rayleigh_angle;
        
        % calculate the location of the rat by both leds
%         pos_average_x = (Cell.pos.x + Cell.pos.x2)/2;
%         pos_average_y = (Cell.pos.y + Cell.pos.y2)/2;
%         
%         % calculate MD
%         [MD_rate_hist,MD_rayleigh_score,MD_rayleigh_angle]=ComputeMovingDirectionality...
%             (Cell.pos.t,pos_average_x,pos_average_y,Cell.spk.t,parms);
%         pop_rayleighMD_score(count) = MD_rayleigh_score;
%         pop_rayleighMD_angle(count) = MD_rayleigh_angle;
%         Cell.MD.rate_hist = MD_rate_hist;
%         Cell.MD.rayleigh.score = MD_rayleigh_score;
%         Cell.MD.rayleigh.angle = MD_rayleigh_angle;
%  
        % insert variance and mean of gc size
        Cell.variance_begin= S_begin.variance;
        Cell.variance_end= S_end.variance;
        
        Cell.mean_values_begin = S_begin.mean_values_list;
        Cell.mean_values_end = S_end.mean_values_list;
        
        Cell.sorted_means_begin=S_begin.sorted_means;
        Cell.sorted_means_end=S_end.sorted_means;
        Cell.sorted_means_total =S_total.sorted_means; 
        
        Cell.norm_means_begin = S_begin.norm_means;
        Cell.norm_means_end = S_end.norm_means;
        
        Cell.peak_rates_begin = S_begin.peak_rates_list;
        Cell.peak_rates_end = S_end.peak_rates_list;
        
        Cell.max_index_begin = S_begin.max_index ;
        Cell.max_index_end  = S_end.max_index; 
        
        Cell.norm_max_index_begin= S_begin.norm_max_index;
        Cell.norm_max_index_end= S_end.norm_max_index;
        
        Cell.zone_mat_begin = S_begin.zone_mat;
        Cell.zone_mat_end=S_end.zone_mat;
        
        Cell.number_zone_mat_begin= S_begin.number_zone_mat;
        Cell.number_zone_mat_end= S_end.number_zone_mat;
        
        Cell.peak_zone_mat_begin= S_begin.peak_zone_mat;
        Cell.peak_zone_mat_end= S_end.peak_zone_mat;
        
        Cell.max_peak_dist_begin= S_begin.max_peak_distance;
        Cell.max_peak_dist_end= S_end.max_peak_distance;
        
        Cell.max_inds_begin = S_begin.max_inds;
        Cell.max_inds_end=S_end.max_inds;
        
        Cell.auto_max_inds_begin= S_begin.auto_max_inds; 
        Cell.auto_max_inds_end= S_end.auto_max_inds; 
        
        Cell.autocorr_begin = S_begin.rate_mat_after_auto_pearson;
        Cell.autocorr_end = S_end.rate_mat_after_auto_pearson;
        
        Cell.six_orientation_pts_begin= S_begin.six_orientation_pts;
        Cell.six_orientation_pts_end= S_end.six_orientation_pts;
        
          Cell.overall_mean_begin= S_begin.overall_mean;
            Cell.overall_mean_end= S_end.overall_mean;
        
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
            m = 4;
            subplot(n,m,1);
            plot(pos_mean_x_begin,pos_mean_y_begin,'k');hold on;
            plot(spk_x_begin,spk_y_begin,'.r');
            axis equal;axis off;
            axis ij
            title(file_name);
            
            subplot(n,m,5);
            plot(pos_mean_x_end,pos_mean_y_end,'k');hold on;
            plot(spk_x_end,spk_y_end,'.r');
            axis equal;axis off;
            axis ij
            
            subplot(n,m,9);
            plot(pos_mean_x,pos_mean_y,'k');hold on;
            plot(spk_x,spk_y,'.r');
            axis equal;axis off;
            axis ij
            
            subplot(n,m,2);
            imagesc(Cell.rate_mat_begin);
            axis image 
            hold on;
    
            subplot(n,m,6);
            imagesc(Cell.rate_mat_end);
            axis image 
            hold on;
    
            subplot(n,m,10);
            imagesc(Cell.rate_mat_total);
            axis image 
            hold on;
            
            subplot(n,m,3)
            imagesc(Cell.zone_mat_begin);
            axis image;
            colorbar;
            
            subplot(n,m,7)
            imagesc(Cell.zone_mat_end);
            axis image;
            colorbar;
            
            subplot(n,m,11)
            imagesc(S_total.zone_mat);
            axis image;
            colorbar;
            
            subplot(n,m,4);
            plot(Cell.sorted_means_begin, 'o-') 
            ymax= max(S_begin.sorted_means);
            xmax= length(S_begin.sorted_means)+1;
            axis ([0 xmax 0 ymax])
            hold on;
            strTitle= sprintf('avg mean=%0.3f \n variability = %0.3f \n mean norm mean=%0.3f \n fano factor= %0.3f', ...
                Cell.overall_mean_begin, Cell.variance_begin, mean(Cell.norm_means_begin), S_begin.var_over_mean);
            title(strTitle);
            
            subplot(n,m,8);
            plot(Cell.sorted_means_end, 'o-') 
            ymax= max(S_end.sorted_means);
            xmax= length(S_end.sorted_means)+1;
            axis ([0 xmax 0 ymax])
            hold on;
            strTitle= sprintf('avg mean=%0.3f \n variability = %0.3f \n mean norm mean=%0.3f \n fano factor= %0.3f', ...
                Cell.overall_mean_end, Cell.variance_end, mean(Cell.norm_means_end), S_end.var_over_mean);
            title(strTitle);
   
           subplot(n,m,12);
            plot(S_total.sorted_means, 'o-') 
            ymax= max(S_total.sorted_means);
            xmax= length(S_total.sorted_means)+1;
            axis ([0 xmax 0 ymax])
            hold on;
            strTitle= sprintf('avg mean=%0.3f \n variability = %0.3f \n mean norm mean=%0.3f \n fano factor= %0.3f',...
                S_total.overall_mean, S_total.variance, mean(S_total.norm_means), S_total.var_over_mean);
            title(strTitle);
            
            
   
cd(parms.dir_save_pictures);
       saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
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
 


function R = get_zones(rate_mat, R)
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
  
%find maximum points in rate_mat_after_auto_pearson 
 
[size_x, size_y] = size(Cell.autocorr);
auto_max_inds_len = 0;
 
%find autocorrelation mat maximum points

for fig_i = 2:size_x-1
    for j = 2:size_y-1
       if Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i+1,j) && ...
           Cell.autocorr(fig_i,j) >Cell.autocorr(fig_i-1,j) && ...
                 Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i,j+1) && ...
                 Cell.autocorr(fig_i,j) > Cell.autocorr(fig_i,j-1)
                   % plot(j,fig_i,'x');
                 
              auto_max_inds_len = auto_max_inds_len+1;    
              auto_max_inds(auto_max_inds_len,1) = fig_i;   %indices of maximum pts
              auto_max_inds(auto_max_inds_len,2) = j;
       end
    end
end


            

% find radius of PF using half distance between center and closest max pt 
 
% cen = 1:auto_max_inds_len;
% distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
%  
% new_distances = sort(distances (:));
% PF_radius = new_distances(2)/2; %% finds 2nd min since min is 0
% 
% % has to be at least 70% of calculated radius distance
% 
% PF_radius = PF_radius;
% 
% S.PF_radius = PF_radius;


% find distances that are too close together between peaks

h=1;
too_close=[];

for cen = 1:max_inds_len-1
   for cen2= (cen+1):max_inds_len
    peak_distance = Distance(max_inds(cen,1), max_inds(cen,2),max_inds(cen2,1), max_inds(cen2,2));
    if peak_distance < 1.25*(Cell.PF_radius) %%%%% add /2 if nonsmooth 
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
if rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) >...    
     rate_mat((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ; 
       
    remove(h) = too_close(cen,2);
 
    h= h+1;
     
elseif rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) <=...    
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
 
[max_inds_len, ~] = size(max_inds);

zone_mat = zeros(size_x, size_y);
 
zone_num = 1;
 
for cen = 1:max_inds_len;
     zone_mat((max_inds(cen,1)),(max_inds(cen,2))) = zone_num;
     
     zone_num = zone_num+1;
end 

peak_zone_mat = zone_mat; 

% if distance to max pt is less than 70% PF radius, assign it value at max pt

R.PF_radius = (0.65* Cell.PF_radius);  %%% ADD /2 if nonsmooth

[size_x, size_y] = size(rate_mat); 

for cen=1:max_inds_len;
     for fig_i =1:size_x
         for j =1:size_y 
             if Distance(fig_i, j, max_inds(cen,1), max_inds(cen,2)) < R.PF_radius        
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

R.max_inds=max_inds;
          
%bkgd_firing = mean_values_list(1);

%S.bkgd_firing = bkgd_firing;

%mean_values_list(1) = []; %drop background mean firing rate

R.number_of_PF = length(peak_rates_list);

R.overall_mean = mean(peak_rates_list); %overall average mean firing rate of cell

R.mean_values_list= peak_rates_list';
%  
R.variance = var(peak_rates_list);

PF_sum= nnz(number_zone_mat);

PF_coverage = PF_sum/(size_x*size_y);

R.PF_coverage = PF_coverage;

R.zone_mat = zone_mat;

% plot mean values
sorted_means = sort(peak_rates_list);

R.sorted_means=sorted_means;

norm_means = R.sorted_means/max(R.sorted_means);

R.norm_means = norm_means;

R.variance_norm_means = var(norm_means);

%storing vars back to S:
 R.new_rate_mat=rate_mat;
 R.size_x=size_x;
 R.size_y=size_y;

 var_over_mean = var(sorted_means)/mean(sorted_means);
R.var_over_mean = var_over_mean;
 
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

R.border_rate_diff = border_to_nonborder_diff;

top_to_central = mean(peak_rates_list(top))/mean(nonborder_means);
R.top = top_to_central;

bottom_to_central = mean(peak_rates_list(bottom))/mean(nonborder_means);
R.bottom = bottom_to_central;

left_to_central= mean(peak_rates_list(left))/mean(nonborder_means);
R.left= left_to_central;

right_to_central= mean(peak_rates_list(right))/mean(nonborder_means);
R.right = right_to_central;

%% finds position of highest firing rate zone

zone_num = find(peak_rates_list==max(peak_rates_list));
zone_num = zone_num(1,1);

location1 = any(top==zone_num);
location2 = any(bottom==zone_num);
location3 = any(left==zone_num);
location4 = any(right==zone_num);

location(1)= location1;
location(2)= location2;
location(3) = location3;
location(4) = location4; 

where = sum(location);

R.location = location; 

R.where = where;

R.peak_rates_list = peak_rates_list;


max_firing_index(1,1) = max_inds(zone_num, 1);
max_firing_index(1,2) = max_inds(zone_num, 2);

R.max_index = max_firing_index;

[size_x, size_y] = size(rate_mat);

norm_max_index(1,1) = max_firing_index(1,1) / size_x;

norm_max_index(1,2) = max_firing_index(1,2) / size_y;

R.norm_max_index = norm_max_index;

R.zone_mat = zone_mat;

R.number_zone_mat = number_zone_mat ;

R.peak_zone_mat = peak_zone_mat;

R.i = i;

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

R.max_peak_distance = max_peak_distance;

disp('');

end %of getzones

            

disp('');
 
end % of main script
 
  
 

 
 
 
 
 
 


