function [norm_shift]= FindPhaseShift(xcorr, spacing)

auto_max_inds= FindAutoMaxInds(xcorr);

xcorr_size=size(xcorr);
center= xcorr_size/2;

dists=nan(1,length(auto_max_inds));
for k=1:length(auto_max_inds)
    dists(k)= Distance(center(1),center(2),auto_max_inds(k,1),...
        auto_max_inds(k,2));
end

phase_shift= min(dists);
norm_shift= phase_shift/mean(spacing);