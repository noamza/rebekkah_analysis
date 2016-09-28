function [ indices, pairs, distances ] = NoamFindPairsOfClosestPoints(a, b)
%PARAMS a and b are sets of point coordinates
%NOAMFINDCLOSESTPEAKS finds set of closest pairs and their distances
indices = []; pairs = []; distances = [];
t = [];
%finds all distances between all combination of points
for i = 1:length(a(:,1)) %for correcting length() function 
    for j = 1:length(b(:,1))%for correcting length() function 
        d = pdist( [a(i,:); b(j,:)] ,'euclidean'); 
        t(end + 1,:) = [d, i, j];
    end
end
%sort point pairs by distance
t = sortrows(t);
%make sure that each point is only used in one pair (by index)
usedA = [-1]; %used from a
usedB = [-1]; %used from b
for i = 1:length(t);
    ia = t(i,2);
    ib = t(i,3);
    %if point hasn't been used, add pair, sorted by distance
    if ~any(ia==usedA) && ~any(ib==usedB) %neither point previously used
        usedA(end + 1) = ia; usedB(end + 1) = ib; 
        distances(end + 1) = t(i,1);
        indices(end + 1,:) = [ia, ib];
        pairs(end +1,:) = [a(ia,:) b(ib,:)];
    end
end
distances = distances'; %why not






