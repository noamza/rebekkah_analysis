function [min_i, min_j] = FindTransfPtToPt2(max_inds1, max_inds2)

[len1, ~]= size(max_inds1);
[len2, ~]= size(max_inds2);

dist=inf;
for i=1:len1
    for j=1:len2
        
        x=max_inds1(i,1);
        y=max_inds1(i,2);
        x1=max_inds2(j,1);
        y1=max_inds2(j,2);
        
        compare_dist= Distance(x,y,x1,y1);
        
        if min(dist, compare_dist) == compare_dist
            min_i=i;
            min_j=j;
            
            dist=compare_dist;
        end
        
    end
end

disp('')