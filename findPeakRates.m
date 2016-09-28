function peak_rates = findPeakRates(max_inds, rate_mat)

try %Noam
    
    [len,~]=size(max_inds);
    peak_rates = nan(1,len);
    for cen = 1:len
        peak_rates(cen)= rate_mat(max_inds(cen,1), max_inds(cen,2));
    end

catch MException %NOAM
    peak_rates = [0, 0];
    warning('findPeakRates: Could not find peak_rates, using dummy values');
end

