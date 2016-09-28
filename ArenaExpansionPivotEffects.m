

%load page image
%extract plot image
%find peaks
%remove peaks that are too close
%find value at peaks
%create figure of 5 images of max peak plotted and peak firing over overall
%...firing ratio
% plot figure of 5 arenas with superimposed corrdinates of normalized max
% ...peak

%load page image

cd('C:\Users\Dori\Dropbox\Dori Rebekkah')
im= imread('barry jeffey NN 2007 grid cell rescaling suppl 2_Page_14.jpg'); %% edit this to how it should be written

%extract plot image

%image 1
[image]= extract_images_from_jpeg(im);
[size_x_1, size_y_1]=size(image);
[max_inds_1, peak_values_1, norm_max_ind1] = findMaxPeaksJPEGImage(image);


%image2
[image_2]= extract_images_from_jpeg(im);
[size_x_2, size_y_2]=size(image_2);
[max_inds_2 , peak_values_2, norm_max_ind2] = findMaxPeaksJPEGImage(image_2);

%image3

[image_3]= extract_images_from_jpeg(im);
[size_x_3, size_y_3]=size(image_3);
[max_inds_3, peak_values_3, norm_max_ind3] = findMaxPeaksJPEGImage(image_3);

%image4
[image_4]= extract_images_from_jpeg(im);
[size_x_4, size_y_4]=size(image_4);
[max_inds_4, peak_values_4, norm_max_ind4] = findMaxPeaksJPEGImage(image_4);

%image5
[image_5]= extract_images_from_jpeg(im);
[size_x_5, size_y_5]=size(image_5);
[max_inds_5, peak_values_5, norm_max_ind5] = findMaxPeaksJPEGImage(image_5);

% figure;
% plot(norm_max_ind1(2), norm_max_ind1(1), 'x', 'MarkerSize', '15'); hold on;
% plot(norm_max_ind2(2), norm_max_ind2(1), 'x', 'MarkerSize', '15'); hold on;
% plot(norm_max_ind3(2), norm_max_ind3(1), 'x', 'MarkerSize', '15'); hold on;
% plot(norm_max_ind4(2), norm_max_ind4(1), 'x', 'MarkerSize', '15'); hold on;
% plot(norm_max_ind5(2), norm_max_ind5(1), 'x', 'MarkerSize', '15'); hold on;
% axis([0 1 0 1]);


fig= figure;
plot(norm_max_ind1(2), norm_max_ind1(1), 'x', 'MarkerSize', 15, 'color', 'r', 'linewidth', 6); hold on;
plot(norm_max_ind2(2), norm_max_ind2(1), 'x', 'MarkerSize', 15, 'linewidth', 6); hold on;
plot(norm_max_ind3(2), norm_max_ind3(1), 'x', 'MarkerSize', 15, 'linewidth', 6); hold on;
plot(norm_max_ind4(2), norm_max_ind4(1), 'x', 'MarkerSize', 15, 'linewidth', 6); hold on;
plot(norm_max_ind5(2), norm_max_ind5(1), 'x', 'MarkerSize', 15, 'color', 'g', 'linewidth', 6); hold on;
axis([0 1 0 1]);
axis ij;

cd('N:\users\rebekkah\barry jeffery pivot point\images')
saveas(fig, 'page14 ex 2.jpg') 

close all 
clearvars