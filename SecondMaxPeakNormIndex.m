function [norm_second_max] = SecondMaxPeakNormIndex(sorted_means, peak_rates, peak_zone_mat)

 second_peak = sorted_means(end-1);
 
  zone_num = find(peak_rates==second_peak);
  
  zone_num = zone_num(1); % if same firing rate at two locations, just take first for analysis
      
 [row, col] = find(peak_zone_mat == zone_num); 

index = [];
index(1,1) = row;
index(1,2) = col;

 [size_x, size_y] = size(peak_zone_mat);
 
  norm_second_max(1,1) = index(1,1) / size_x;
        
  norm_second_max(1,2) = index(1,2) / size_y;
  
  