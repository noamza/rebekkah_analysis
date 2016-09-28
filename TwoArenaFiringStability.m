function [corr_two]=TwoArenaFiringStability(rate_map1,rate_map2,PF_rad1,PF_rad2)

max_inds1=FindMaxIndsRateMap(rate_map1);
peak_rates1 = findPeakRates(max_inds1, rate_map1);
[zm_1, ~]= CreateZoneMat(size(rate_map1), PF_rad1, max_inds1, peak_rates1);

max_inds2=FindMaxIndsRateMap(rate_map2);

if exist('max_inds2','var')
    
    peak_rates2 = findPeakRates(max_inds2, rate_map2);
    [zm_2, ~]= CreateZoneMat(size(rate_map2), PF_rad2, max_inds2, peak_rates2);
    
    [~,zm_2]=ArenaSameSize(zm_1,zm_2);
    
    if length(max_inds1) >= 3
        
        % find peak rates of second session using 1st session inds
        peak_rates_2=findPeakRates(max_inds1,zm_2);
        peak_rates_1=peak_rates1;
        
        prb_orig=peak_rates_1;
        pre_orig=peak_rates_2;
        
        peak_rates_1(prb_orig==0 | pre_orig==0)= [];
        peak_rates_2(prb_orig==0 | pre_orig==0)= [];
        
        if length(peak_rates_1) >2
            corr_two= corrcoef(peak_rates_1,peak_rates_2);
            corr_two=corr_two(2);
        else
            corr_two=nan;
        end
        
        %                     figure;
        %                     subplot(1,2,1)
        %                     imagesc(zm_1)
        %                     subplot(1,2,2)
        %                     imagesc(zm_2)
        
    else
        corr_two=nan;
        
    end        
else corr_two =nan;
        close all
  
end

