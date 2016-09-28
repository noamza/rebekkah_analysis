function [size_peaks]= FindSizeofField(peak_rates, number_zone_mat, peak_zone_mat)


size_peaks= nan(1, length(peak_rates));
    for h=1:length(peak_rates)
        size_peaks(h)=sum(number_zone_mat(number_zone_mat==h)>= 0.9*peak_zone_mat(h));
    end