function extract_images_from_jpeg(im)

min_image_size = 10; % min. number of pixels in an image

figure;
image(im);

white_matrix = im(:,:,1) > 250 & im(:,:,2) > 250 & im(:,:,3) >250;
black_matrix = ~white_matrix;


figure;
imagesc(black_matrix);
    
cc = bwconncomp(black_matrix)

bw = bwselect(black_matrix);

inds  = find(bw);

[x,y] = ind2sub(size(bw),inds);

x_ind = x - min(x) + 1;
y_ind = y - min(y) + 1;


for i = 1:length(x)   
    new_mat(x_ind(i),y_ind(i),:) = im(x(i),y(i),:);
end

black_inds = find(new_mat(:,:,1)  > 200 & new_mat(:,:,2) > 200 & new_mat(:,:,3) > 200);

for i = 1:3
    tmp_im = new_mat(:,:,i);
    tmp_im(black_inds) = 0;
    new_mat(:,:,i) = tmp_im;
end

[n,m,rgb] = size(new_mat);

% get rid of all-zero rows and columns
row_list = [];
for i = 1:n  
   row = new_mat(i,:,1);
   if ~all(row == 0)
       row_list(end+1) = i;
   end
end

 col_list = [];
for j = 1:m 
   col = new_mat(:,j,1);
   if ~all(col == 0)
       col_list(end+1) = j;
   end
end  

cropped_mat = new_mat(row_list,col_list,:);



% find the closest values in the 'jet' colormap

cmap = colormap(jet(256));
cmap = cmap*256;

mat = cropped_mat;

for i = 1:size(mat,1);
    for j = 1:size(mat,2);
        
        rgb_values = squeeze(mat(i,j,:));
        
        comp_vec = double(cmap) - repmat(double(rgb_values'),[256 1]);
        [val,ind] = min(sum(comp_vec.*comp_vec,2));
        out_mat(i,j) = ind;
    end
end

        
        
        


