function S= get_zones(rate_mat, S, PF_radius)
        %find max points in rate mat
        %find max points in auto mat, use this to find radius of PF
        %use PF radius in rate mat to create zone mat by comparing max pts to PF
        %radius
        
        %% works same as smooth get_zones except need to input PF radius size to work (unable to create it itself due to noise of autocorr)
        
        
        rate_mat_orig = rate_mat;
        rate_mat(isnan(rate_mat))=0;    %change nans to zeros
        
        [size_x,size_y]=size(rate_mat);
        
        max_inds = FindMaxIndsRateMap(rate_mat); % look for local maxima in rate_mat
        
        rate_mat = rate_mat_orig;       %return orignal rate mat without zeros border
        [size_x,size_y]=size(rate_mat);
        
        %create auto correlated rate map
        
        % autocorr = Cross_Correlation(rate_mat,rate_mat);
     
        % S.autocorr= autocorr;
        
        
        %find maximum points in rate_mat_after_auto_pearson
        
     %   auto_max_inds= FindAutoMaxInds(autocorr);
        
      %  S.auto_max_ins= auto_max_inds;
        % find radius of PF using half distance between center and closest max pt
        
      % PF_radius= findPlaceFieldRadius(autocorr, auto_max_inds);
      % S.PF_radius= PF_radius; 
        
      
    %   auto_max_inds= RemoveTooCloseMaxInds(auto_max_inds, PF_radius, autocorr);
             
     %   PF_radius= findPlaceFieldRadius(autocorr, auto_max_inds);
     %   S.PF_radius= PF_radius;         

        % find peaks that are too close to eachother
        max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, 1.8);
        
   
        % give each max point a value
        
        [size_x, size_y] = size(rate_mat);
        
        [max_inds_len, ~]= size(max_inds);
        
        zone_mat = zeros(size_x, size_y);
        
        zone_num = 1;
        
        for cen = 1:max_inds_len;
            zone_mat((max_inds(cen,1)),(max_inds(cen,2))) = zone_num;
            
            zone_num = zone_num+1;
        end
        
        
        peak_zone_mat = zone_mat;
        
        % if distance to max pt is less than 70% PF radius, assign it value at max pt
        
        [size_x, size_y] = size(rate_mat);
        
        for cen=1:max_inds_len;
            for fig_i =1:size_x
                for j =1:size_y
                    if Distance(fig_i, j, max_inds(cen,1), max_inds(cen,2)) < 0.8 * PF_radius  %change this depending on how large you want fields to be
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
        
%         rate_mat(isnan(rate_mat))= 0;
%         
%         peak_rates_len = 0;
%         
%         for cen = 1:number_zones_len;
%         
%               find_spk = find(number_zone_mat == number_zones(cen));
%               if sum(find_spk>length(rate_mat))>0
%                 disp('')
%               end
%         
%               peak_rates_len= peak_rates_len+1;
%               peak_rates_list(peak_rates_len) = mean(rate_mat(find_spk));
%         end
        
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
        
%        
%         
        %% finds position of highest firing rate zone
        
         zone_num = find(peak_rates_list==max(peak_rates_list));
%         
        zone_num= zone_num(1) ;     %if there are multiple max fields of same firing rate, take the 1st

%         location(1) = any(left==zone_num);
%         location(2) = any(right==zone_num);
%         location(3) = any(bottom==zone_num);
%         location(4) = any(top==zone_num);
%         
%         where = sum(location);
%         
%         S.location = location;
%         
%         S.where = where;
        
        S.peak_rates = peak_rates_list;
        
        
        max_firing_index(1) = max_inds(zone_num, 1);
        max_firing_index(2) = max_inds(zone_num, 2);
        
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
        
     %   [size_x, size_y] = size(S.autocorr);
        
    %    auto_distances = [];
    %    new_distances = [];
        
        %%finds distances from center of auto_corr_map to all max peaks
       
        
 %       cen = 1:length(auto_max_inds);
  %      auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
        
%        new_distances = sort(auto_distances (:));
        
        
        % is the max peak stronger firing than the 2nd max peak by at least
        % 90% of its firing
        
%         if sorted_means(end-1) > sorted_means(end)*0.87
%             pivot= 'not strong'; 
%         elseif sorted_means(end-1) < sorted_means(end)*0.87
%             pivot= 'strong';
%         end
% 
%         S.pivot=pivot;
        
      disp('');  
        
    end %of getzones