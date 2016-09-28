function [phase_dist]= FindPhaseDiff(acorr_xcorr_hex, xcorr)

center_pt= acorr_xcorr_hex(7,:);

[max_peaks] = FindAutoMaxInds(xcorr);

dists_center=nan(1,length(max_peaks));
for h= 1:length(max_peaks)
dists_center(h)= Distance(max_peaks(h,1), max_peaks(h,2),...
                            center_pt(1), center_pt(2));
end
                        
min_dist= min(dists_center);
ind= find(dists_center(dists_center==min_dist));
ind=ind(1);

phase_dist= dists_center(ind); 