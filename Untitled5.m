function PrintImages

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat')

for i=1:length(zone_mats_all)

fig=figure;
n=2;m=3;
subplot(n,m,1)
plot(pos_x_all{i},pos_y_all{i},'k');hold on;
plot(spk_x_all{i},spk_y_all{i},'.r');
axis equal;axis off;
axis ij;
title(sprintf('cell %d', i));

subplot(n,m,2)
imagesc(rate_mats_all{i});
axis equal; axis off;
title(sprintf('%0.1f Hz', max(peak_rates_sm_all{i})), 'HorizontalAlignment', 'left');

subplot(n,m,3)
imagesc(zone_mats_all{i});
axis equal; axis off;

subplot(n,m,4)
plot(1:length(peak_rates_sm_all{i}),sort(peak_rates_sm_all{i}), 'ko-');
title(sprintf('Fano factor=%0.1f', fanos(i)));

subplot(n,m,5)
imagesc(autocorrs_all{i});
axis equal; axis off;

subplot(n,m,6)
fourier_mat=fft2(rate_mats_all{i});
fourier_mat(1)=0;
[x,y] = size(fourier_mat);
fourier_mat = fftshift(fourier_mat);
x1 = round(x/4); x2 = round(3*x/4); y1 = round(y/4); y2 = round(3*y/4);
fourier_mat = fourier_mat(x1:x2,y1:y2);
imagesc(abs(fourier_mat)); 
axis equal;axis off;

%save images
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7\images');
saveas(fig,sprintf('%d.jpg',i));


close all
end