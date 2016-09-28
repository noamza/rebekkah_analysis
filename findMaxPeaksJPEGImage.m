function [max_inds, value_list, norm_max_ind] = findMaxPeaksJPEGImage(image_1) 

% look for local maxima in image
 
max_inds_len = 0;
max_inds = [];

[size_x size_y] = size(image_1);



image_rate= nan([size_x+2, size_y+2]); 
image_rate(2:size_x+1, 2:size_y+1)=image_1(1:size_x, 1:size_y);



for i = 2:size_x-1
    for j = 2:size_y-1       
        if image_rate(i,j) > image_rate(i+1,j) && ...
                image_rate(i,j) > image_rate(i-1,j) && ...
                image_rate(i,j) > image_rate(i,j+1) && ...
                image_rate(i,j) > image_rate(i,j-1)
                hold on
            max_inds_len = max_inds_len+1;    
            max_inds(max_inds_len,1) = i-1;     %list of maximum pt indices
            max_inds(max_inds_len,2) = j-1;
        end
    end
end


% find distances that are too close together between peaks

h=1;
too_close=[];

for cen = 1:max_inds_len-1
   for cen2= (cen+1):max_inds_len
    peak_distance = Distance(max_inds(cen,1), max_inds(cen,2),max_inds(cen2,1), max_inds(cen2,2));
    if peak_distance < 20
     too_close (h,1) = cen;
     too_close(h,2) = cen2;
    
     h=h+1;       
    end
   end
end

if ~isempty (too_close)

    [size_x, size_y] = size(too_close);
    
    too_close_len = size_x;


%remove one of the peaks that are too close with lower firing rate

h=1;

for cen = 1: too_close_len
if image_rate((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) >...    
     image_rate((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ; 
       
    remove(h) = too_close(cen,2);
 
    h= h+1;
     
elseif image_rate((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) <=...    
     image_rate((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ;

    remove(h)= too_close(cen,1);

    h=h+1;
    
else
    disp ('wtf')
end
end

    max_inds(remove,1) = 0;
    max_inds(remove,2) = 0;
    
max_inds(all(max_inds==0,2),:) = [];
 
end 

for cen=1:length(max_inds)
    value_list(cen)= image_1(max_inds(cen,1), max_inds(cen,2));
end

[max_value, max_ind]= max(value_list);

[size_x size_y] = size(image_1);
norm_max_ind(1)= max_inds(max_ind,1)/ size_x;
norm_max_ind(2)= max_inds(max_ind,2)/ size_y;

end 