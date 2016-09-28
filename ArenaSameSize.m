function [zone_mat_1, zone_mat_2]= ...
    ArenaSameSize(zone_mat_1, zone_mat_2)

size_1= size(zone_mat_1);
size_2= size(zone_mat_2);

new_size= [max(size_1(1),size_2(1)) max(size_1(2),size_2(2))];

% stretch so same size (rate mats a bit off)
if ~isequal(size_1, new_size)
    [zone_mat_1] = StretchImage(zone_mat_1, size_1, new_size);
end

if ~isequal(size_2, new_size)
    [zone_mat_2] = StretchImage(zone_mat_2, size_2, new_size);
end

