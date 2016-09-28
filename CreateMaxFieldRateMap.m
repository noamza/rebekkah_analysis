function [max_pt_rate_map]= CreateMaxFieldRateMap(total_max_inds)

max_pt_rate_map= zeros(20,20);

total_max_inds_conv= round((total_max_inds/5) * 100);

total_max_inds_conv(total_max_inds_conv==0) = 1;

for h = 1:length(total_max_inds)
        max_pt_rate_map(total_max_inds_conv(h,1),total_max_inds_conv(h,2))= max_pt_rate_map(total_max_inds_conv(h,1),total_max_inds_conv(h,2))+1;
end
