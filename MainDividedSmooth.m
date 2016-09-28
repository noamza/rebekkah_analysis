function Main
% Date: 24 of July 2014
dbstop if error

parms.dir_load_data = 'C:\Users\Dori\Desktop\elipse plots PF 8 results';

parms.dir_save_data = 'N:\users\rebekkah\does fanofactor decrease over time\results 3 4 and 6 timestamps';

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
    
    if pos_t_len >= 45000 & pos_t_len < 60000
        
        one_inx= 1:15000;
        two_inx= 15001:30000;
        three_inx= 30001:45000;
        
        
        pos_t_1= Cell.pos.t(one_inx);   %first half of session
        pos_t_2= Cell.pos.t(two_inx); %second half of session
        pos_t_3= Cell.pos.t(three_inx);
        
        tmp_1 = abs(Cell.spk.t-pos_t_1(end)); %subtract last timestamp of first session
        tmp_2= abs(Cell.spk.t- pos_t_2(end)) ;% subtract first timestamp of second session
        tmp_3= abs(Cell.spk.t- pos_t_3(end)) ;
        [~, idx_1] = min(tmp_1); %index of first session
        [~, idx_2] = min(tmp_2) ;%index of first session
        [~, idx_3] = min(tmp_3);
        spk_t_1 = Cell.spk.t(1:idx_1); %spike vector of first sesh
        spk_t_2 = Cell.spk.t(idx_1+1:idx_2); %spike vector of second sesh
        spk_t_3 = Cell.spk.t(idx_2+1:idx_3) ;
        
        if ~isempty(Cell.pos.x2)
            dt = Cell.pos.t(2)-Cell.pos.t(1);
            
            % calculate the the rat's head direction (using average of both leds)
            pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
            pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
            
            % build the axis of location when spikes are made
            spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
            spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
            
            % divide into first 5 minutes and last 5 minutes
            pos_mean_x_1= pos_mean_x(one_inx); %x position for first half of session
            pos_mean_x_2= pos_mean_x(two_inx); %x position for second half of session
            pos_mean_x_3= pos_mean_x(three_inx);
            pos_mean_y_1= pos_mean_y(one_inx); %y position for first half of session
            pos_mean_y_2= pos_mean_y(two_inx); %y position for second half of session
            pos_mean_y_3= pos_mean_y(three_inx);
            
            % confirm the Cell.spk.t and spk_x are same size and correlate to
            % eachother!!!!!!!!!!!!!!!
            
            spk_x_1= spk_x(1:idx_1);
            spk_x_2= spk_x(idx_1+1:idx_2);
            spk_x_3= spk_x(idx_2+1:idx_3);
            spk_y_1= spk_y(1:idx_1);
            spk_y_2= spk_y(idx_1+1:idx_2);
            spk_y_3= spk_y(idx_2+1:idx_3);
            
            % get rate matrix - using function: Creat_Rate_Map
            rate_mat_1=CreateRateMap(pos_mean_x_1,pos_mean_y_1,pos_t_1,spk_x_1,spk_y_1,spk_t_1,parms);
            rate_mat_2=CreateRateMap(pos_mean_x_2,pos_mean_y_2,pos_t_2,spk_x_2,spk_y_2,spk_t_2,parms);
            rate_mat_3=CreateRateMap(pos_mean_x_3,pos_mean_y_3,pos_t_3,spk_x_3,spk_y_3,spk_t_3,parms);
            
            S_1=get_zones(rate_mat_1, dat.S);
            S_2=get_zones(rate_mat_2, dat.S);
            S_3=get_zones(rate_mat_3, dat.S);
            
            fanofactor_all=[];
            fanofactor_all(1)= var(S_1.sorted_means)/mean(S_1.sorted_means);
            fanofactor_all(2)= var(S_2.sorted_means)/mean(S_2.sorted_means);
            fanofactor_all(3)= var(S_3.sorted_means)/mean(S_3.sorted_means);
            
            
            timestamps= [1 2 3];
            
            corr_fanofactor= corrcoef(fanofactor_all, timestamps);
            
            Cell.corr_fanofactor= corr_fanofactor;
            
            corr_fanofactor_all(i)= corr_fanofactor(2);
            
            %  Cell_num(count)=i;
            
            
            
            %  count=count+1;
            % alpha
            %      close all;
        end
        
        % saves results
        cd(parms.dir_save_data);
        title1 =strcat (file_name);
        %title1 =strcat ('Results_',file_name);
        save (title1,'S_1', 'S_2', 'S_3');
        
    elseif pos_t_len >= 60000  & pos_t_len < 75000
        
        one_inx= 1:15000;
        two_inx= 15001:30000;
        three_inx= 30001:45000;
        four_inx= 45001:60000;
        
        
        pos_t_1= Cell.pos.t(one_inx);   %first half of session
        pos_t_2= Cell.pos.t(two_inx); %second half of session
        pos_t_3= Cell.pos.t(three_inx);
        pos_t_4= Cell.pos.t(four_inx);
        
        tmp_1 = abs(Cell.spk.t-pos_t_1(end)); %subtract last timestamp of first session
        tmp_2= abs(Cell.spk.t- pos_t_2(end));% subtract first timestamp of second session
        tmp_3= abs(Cell.spk.t- pos_t_3(end));
        tmp_4= abs(Cell.spk.t- pos_t_4(end));
        [~, idx_1] = min(tmp_1); %index of first session
        [~, idx_2] = min(tmp_2) ;%index of first session
        [~, idx_3] = min(tmp_3);
        [~, idx_4] = min(tmp_4);
        spk_t_1 = Cell.spk.t(1:idx_1); %spike vector of first sesh
        spk_t_2 = Cell.spk.t(idx_1+1:idx_2); %spike vector of second sesh
        spk_t_3 = Cell.spk.t(idx_2+1:idx_3) ;
        spk_t_4 = Cell.spk.t(idx_3+1:idx_4) ;
        
        if ~isempty(Cell.pos.x2)
            dt = Cell.pos.t(2)-Cell.pos.t(1);
            
            % calculate the the rat's head direction (using average of both leds)
            pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
            pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
            
            % build the axis of location when spikes are made
            spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
            spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
            
            % divide into first 5 minutes and last 5 minutes
            pos_mean_x_1= pos_mean_x(one_inx); %x position for first half of session
            pos_mean_x_2= pos_mean_x(two_inx); %x position for second half of session
            pos_mean_x_3= pos_mean_x(three_inx);
            pos_mean_x_4= pos_mean_x(four_inx);
            pos_mean_y_1= pos_mean_y(one_inx); %y position for first half of session
            pos_mean_y_2= pos_mean_y(two_inx); %y position for second half of session
            pos_mean_y_3= pos_mean_y(three_inx);
            pos_mean_y_4= pos_mean_y(four_inx);
            
            % confirm the Cell.spk.t and spk_x are same size and correlate to
            % eachother!!!!!!!!!!!!!!!
            
            spk_x_1= spk_x(1:idx_1);
            spk_x_2= spk_x(idx_1+1:idx_2);
            spk_x_3= spk_x(idx_2+1:idx_3);
            spk_x_4= spk_x(idx_3+1:idx_4);
            spk_y_1= spk_y(1:idx_1);
            spk_y_2= spk_y(idx_1+1:idx_2);
            spk_y_3= spk_y(idx_2+1:idx_3);
            spk_y_4= spk_y(idx_3+1:idx_4);
            
            % get rate matrix - using function: Creat_Rate_Map
            rate_mat_1=CreateRateMap(pos_mean_x_1,pos_mean_y_1,pos_t_1,spk_x_1,spk_y_1,spk_t_1,parms);
            rate_mat_2=CreateRateMap(pos_mean_x_2,pos_mean_y_2,pos_t_2,spk_x_2,spk_y_2,spk_t_2,parms);
            rate_mat_3=CreateRateMap(pos_mean_x_3,pos_mean_y_3,pos_t_3,spk_x_3,spk_y_3,spk_t_3,parms);
            rate_mat_4=CreateRateMap(pos_mean_x_4,pos_mean_y_4,pos_t_4,spk_x_4,spk_y_4,spk_t_4,parms);
            
            S_1=get_zones(rate_mat_1, dat.S);
            S_2=get_zones(rate_mat_2, dat.S);
            S_3=get_zones(rate_mat_3, dat.S);
            S_4=get_zones(rate_mat_4, dat.S);
            
            fanofactor_all(1)= var(S_1.sorted_means)/mean(S_1.sorted_means);
            fanofactor_all(2)= var(S_2.sorted_means)/mean(S_2.sorted_means);
            fanofactor_all(3)= var(S_3.sorted_means)/mean(S_3.sorted_means);
            fanofactor_all(4)= var(S_4.sorted_means)/mean(S_4.sorted_means);
            
            timestamps= [1 2 3 4];
            
            corr_fanofactor= corrcoef(fanofactor_all, timestamps);
            
            Cell.corr_fanofactor= corr_fanofactor;
            
            corr_fanofactor_all(i)= corr_fanofactor(2);
            
            %  Cell_num(count)=i;
            
            
            
            %  count=count+1;
            % alpha
            %      close all;
        end
        
        % saves results
        cd(parms.dir_save_data);
        title1 =strcat (file_name);
        %title1 =strcat ('Results_',file_name);
        save (title1,'S_1', 'S_2', 'S_3', 'S_4');
        
        
    elseif pos_t_len >= 90000 % & pos_t_len < 120000
        
        one_inx= 1:15000;
        two_inx= 15001:30000;
        three_inx= 30001:45000;
        four_inx= 45001:60000;
        five_inx= 60001:75000;
        six_inx= 75001:90000;
        
        
        pos_t_1= Cell.pos.t(one_inx);   %first half of session
        pos_t_2= Cell.pos.t(two_inx); %second half of session
        pos_t_3= Cell.pos.t(three_inx);
        pos_t_4= Cell.pos.t(four_inx);
        pos_t_5= Cell.pos.t(five_inx);
        pos_t_6= Cell.pos.t(six_inx);
        
        tmp_1 = abs(Cell.spk.t-pos_t_1(end)); %subtract last timestamp of first session
        tmp_2= abs(Cell.spk.t- pos_t_2(end));% subtract first timestamp of second session
        tmp_3= abs(Cell.spk.t- pos_t_3(end));
        tmp_4= abs(Cell.spk.t- pos_t_4(end));
        tmp_5= abs(Cell.spk.t- pos_t_5(end));
        tmp_6= abs(Cell.spk.t- pos_t_6(end));
        [~, idx_1] = min(tmp_1); %index of first session
        [~, idx_2] = min(tmp_2) ;%index of first session
        [~, idx_3] = min(tmp_3);
        [~, idx_4] = min(tmp_4);
        [~, idx_5] = min(tmp_5);
        [~, idx_6] = min(tmp_6);
        spk_t_1 = Cell.spk.t(1:idx_1); %spike vector of first sesh
        spk_t_2 = Cell.spk.t(idx_1+1:idx_2); %spike vector of second sesh
        spk_t_3 = Cell.spk.t(idx_2+1:idx_3) ;
        spk_t_4 = Cell.spk.t(idx_3+1:idx_4) ;
        spk_t_5 = Cell.spk.t(idx_4+1:idx_5) ;
        spk_t_6 = Cell.spk.t(idx_5+1:idx_6) ;
        
        if ~isempty(Cell.pos.x2)
            dt = Cell.pos.t(2)-Cell.pos.t(1);
            
            % calculate the the rat's head direction (using average of both leds)
            pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
            pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
            
            % build the axis of location when spikes are made
            spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
            spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
            
            % divide into first 5 minutes and last 5 minutes
            pos_mean_x_1= pos_mean_x(one_inx); %x position for first half of session
            pos_mean_x_2= pos_mean_x(two_inx); %x position for second half of session
            pos_mean_x_3= pos_mean_x(three_inx);
            pos_mean_x_4= pos_mean_x(four_inx);
            pos_mean_x_5= pos_mean_x(five_inx);
            pos_mean_x_6= pos_mean_x(six_inx);
            pos_mean_y_1= pos_mean_y(one_inx); %y position for first half of session
            pos_mean_y_2= pos_mean_y(two_inx); %y position for second half of session
            pos_mean_y_3= pos_mean_y(three_inx);
            pos_mean_y_4= pos_mean_y(four_inx);
            pos_mean_y_5= pos_mean_y(five_inx);
            pos_mean_y_6= pos_mean_y(six_inx);
            
            % confirm the Cell.spk.t and spk_x are same size and correlate to
            % eachother!!!!!!!!!!!!!!!
            
            spk_x_1= spk_x(1:idx_1);
            spk_x_2= spk_x(idx_1+1:idx_2);
            spk_x_3= spk_x(idx_2+1:idx_3);
            spk_x_4= spk_x(idx_3+1:idx_4);
            spk_x_5= spk_x(idx_4+1:idx_5);
            spk_x_6= spk_x(idx_5+1:idx_6);
            spk_y_1= spk_y(1:idx_1);
            spk_y_2= spk_y(idx_1+1:idx_2);
            spk_y_3= spk_y(idx_2+1:idx_3);
            spk_y_4= spk_y(idx_3+1:idx_4);
            spk_y_5= spk_y(idx_4+1:idx_5);
            spk_y_6= spk_y(idx_5+1:idx_6);
            
            % get rate matrix - using function: Creat_Rate_Map
            rate_mat_1=CreateRateMap(pos_mean_x_1,pos_mean_y_1,pos_t_1,spk_x_1,spk_y_1,spk_t_1,parms);
            rate_mat_2=CreateRateMap(pos_mean_x_2,pos_mean_y_2,pos_t_2,spk_x_2,spk_y_2,spk_t_2,parms);
            rate_mat_3=CreateRateMap(pos_mean_x_3,pos_mean_y_3,pos_t_3,spk_x_3,spk_y_3,spk_t_3,parms);
            rate_mat_4=CreateRateMap(pos_mean_x_4,pos_mean_y_4,pos_t_4,spk_x_4,spk_y_4,spk_t_4,parms);
            rate_mat_5=CreateRateMap(pos_mean_x_5,pos_mean_y_5,pos_t_5,spk_x_5,spk_y_5,spk_t_5,parms);
            rate_mat_6=CreateRateMap(pos_mean_x_6,pos_mean_y_6,pos_t_6,spk_x_6,spk_y_6,spk_t_6,parms);
            
            S_1=get_zones(rate_mat_1, dat.S);
            S_2=get_zones(rate_mat_2, dat.S);
            S_3=get_zones(rate_mat_3, dat.S);
            S_4=get_zones(rate_mat_4, dat.S);
            S_5=get_zones(rate_mat_5,dat.S);
            S_6=get_zones(rate_mat_6, dat.S);
            
            
            fanofactor_all(1)= var(S_1.sorted_means)/mean(S_1.sorted_means);
            fanofactor_all(2)= var(S_2.sorted_means)/mean(S_2.sorted_means);
            fanofactor_all(3)= var(S_3.sorted_means)/mean(S_3.sorted_means);
            fanofactor_all(4)= var(S_4.sorted_means)/mean(S_4.sorted_means);
            fanofactor_all(5)= var(S_5.sorted_means)/mean(S_5.sorted_means);
            fanofactor_all(6)= var(S_6.sorted_means)/mean(S_6.sorted_means);
            
            timestamps= [1 2 3 4 5 6];
            
            corr_fanofactor= corrcoef(fanofactor_all, timestamps);
            
            Cell.corr_fanofactor= corr_fanofactor;
            
            corr_fanofactor_all(i)= corr_fanofactor(2);
            
            %  Cell_num(count)=i;
            
            
            
            %  count=count+1;
            % alpha
            %      close all;
        end
        
        % saves results
        cd(parms.dir_save_data);
        title1 =strcat (file_name);
        save (title1,'S_1', 'S_2', 'S_3', 'S_4', 'S_5');
        
        
    elseif pos_t_len < 45000
        
        corr_fanofactor_all(i)= nan;
        
    end
end



%count = count - 1;

save('corrcoeff fanofactor', 'corr_fanofactor_all')

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














