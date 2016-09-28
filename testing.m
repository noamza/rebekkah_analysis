
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info')
load('rescaling arenas data info.mat')

rate_mats=rate_mats_all{4};

figure;
n=5;
m=3;

inc=0.1725;

subplot(n,m,1)
imagesc(rate_mats{1})
axis off;axis square;
annotation('arrow', [0.4 0.6],[0.86 0.86]);

subplot(n,m,4)
imagesc(rate_mats{2})
axis off;axis image;
annotation('arrow', [0.4 0.6],[0.86-inc 0.86-inc]);

subplot(n,m,7)
imagesc(rate_mats{3})
axis off;axis equal;
annotation('arrow', [0.4 0.6],[0.86-(inc*2) 0.86-(inc*2)]);

subplot(n,m,10)
imagesc(rate_mats{4})
axis off;axis image;
annotation('arrow', [0.4 0.6],[0.86-(inc*3) 0.86-(inc*3)]);

subplot(n,m,13)
imagesc(rate_mats{5})
axis off;axis square;
annotation('arrow', [0.4 0.6],[0.86-(inc*4) 0.86-(inc*4)]);

%.....

subplot(n,m,3)
imagesc(rate_mats{1})
axis off;axis square;

subplot(n,m,6)
imagesc(rate_mats{2})
axis off;axis square;

subplot(n,m,9)
imagesc(rate_mats{3})
axis off;axis square;

subplot(n,m,12)
imagesc(rate_mats{4})
axis off;axis square;

subplot(n,m,15)
imagesc(rate_mats{5})
axis off;axis square;

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 30.45], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 30.45])


