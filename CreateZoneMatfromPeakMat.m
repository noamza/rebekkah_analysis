function [number_zone_mat] = CreateZoneMatfromPeakMat(peak_zone_mat, max_inds, PF_radius); 


[size_x, size_y] = size(peak_zone_mat);
        
        for cen=1:length(max_inds);
            for fig_i =1:size_x
                for j =1:size_y
                    if Distance(fig_i, j, max_inds(cen,1), max_inds(cen,2)) < PF_radius  %change this depending on how large you want fields to be
                        number_zone_mat(fig_i,j)= cen;
                    end
                end
            end
        end