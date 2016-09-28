function [rates_corr, orig, next]= FindFiringProfileGC(max_inds_orig, max_inds_next, ...
    rates_orig, rates_next)

% takes two sets of max_indices
% finds closest sets of points
% compares peak_rates from one to the other
% correlation is how stable the firing gradient is between the two
peak_rates_orig=rates_orig;

[len1,~]= size(max_inds_orig);
[len2,~]= size(max_inds_next);

%len_diff= len1/len2;    %ratio diff of PF nums

if len2~= 1 && len1 ~=1 % && abs(1- (1/len_diff)) < 0.5  
    
    min_len= min(len1,len2); % number of least fields
    
    zone_val= nan([1 len1]);
    
    % make sure both hyperfields are included in the corr
         % analysis:BELOW
    % for case with more fields, make sure hyperfield included first         
    if len1 > len2
        hyperind= find(rates_orig==max(rates_orig)); 

        [min_i, min_j] = FindTransfPtToPt2(max_inds_orig(hyperind,:), ...
            max_inds_next);
        
        min_i=hyperind;

        max_inds_orig(min_i,:)= nan;
        max_inds_next(min_j,:)= nan;
        
        zone_val(min_j)= min_i;
    
        min_len=min_len-1;

    elseif len2>len1

        hyperind= find(rates_next==max(rates_next)); 

        [min_i, min_j] = FindTransfPtToPt2(max_inds_orig, ...
            max_inds_next(hyperind,:));
        
        min_j=hyperind;

        max_inds_orig(min_i,:)= nan;
        max_inds_next(min_j,:)= nan;
        
        zone_val(min_j)= min_i;
    
        min_len=min_len-1;
    end
    % make sure both hyperfields are included in corr analysis:ABOVE

    %find pairs of field centers closest together
    for h=1:min_len
        
        [min_i, min_j] = FindTransfPtToPt2(max_inds_orig, max_inds_next);
        
        max_inds_orig(min_i,:)= nan;
        max_inds_next(min_j,:)= nan;
        
        zone_val(min_j)= min_i;
    
      
    end
    
    zone_val(zone_val==0)= nan; % produced a 0 as an ind, added to prevent
    
    peak_rates_next=nan(1,len1);
    for h=1:length(zone_val);
        if ~isnan(zone_val(h))
            peak_rates_next(zone_val(h))=rates_next(h);
        end
    end
    
    zone_missing= setdiff(1:len1, zone_val);
    
    try
        peak_rates_orig(zone_missing)=nan;
    end
    
    nan_ind= find(isnan(peak_rates_orig));
    peak_rates_orig(nan_ind)=[];
    peak_rates_next(nan_ind)=[];
    
    rates_corr= corrcoef(peak_rates_orig, peak_rates_next);
    rates_corr=rates_corr(2);
    
    orig= peak_rates_orig;
    next= peak_rates_next;

else
    rates_corr= nan;
    
    
end




