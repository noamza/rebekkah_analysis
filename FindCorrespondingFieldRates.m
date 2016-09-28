function [peak_rates_next]= FindCorrespondingFieldRates(max_inds_orig,max_inds_next,rates_next,min_len)

len1=length(max_inds_orig);

zone_val= nan([1 len1]);

%find pairs of field centers closest together
for h=1:min_len
    
    [min_i, min_j] = FindTransfPtToPt2(max_inds_orig, max_inds_next);
    
    max_inds_orig(min_i,:)= nan;
    max_inds_next(min_j,:)= nan;
    
    zone_val(min_j)= min_i;
end

%organizes peak rates to corresponding field index
peak_rates_next=nan(1,len1);
for h=1:length(zone_val);
    if ~isnan(zone_val(h))
        peak_rates_next(zone_val(h))=rates_next(h);
    end
end
