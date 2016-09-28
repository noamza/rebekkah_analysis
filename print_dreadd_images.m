figure; subplot(3,4,1)
plot(pos_mean_x,pos_mean_y,'k'); axis image; hold on;
plot(spkx,spky,'.r');
axis off;
title(sprintf('Date: 22.01 \n T3C1 \n Pre-inj. of CNO' ), 'fontname', 'calibri')

subplot(3,4,5)
imagesc(rate_mat); axis xy; axis square; axis off;
title(sprintf('%0.2f', peak_rate), 'fontname', 'calibri', 'HorizontalAlignment', 'left')

subplot(3,4,9)
imagesc(rate_mat); axis xy; axis square; axis off;
title(sprintf('%0.2f', peak_rate), 'fontname', 'calibri', 'HorizontalAlignment', 'left');

lim= peak_rate;

%% 2nd column

subplot(3,4,2)
plot(pos_mean_x,pos_mean_y,'k'); axis image; hold on;
plot(spkx,spky,'.r');
axis off;
title('90 mins post-inj.', 'fontname', 'calibri')

subplot(3,4,6)
imagesc(rate_mat); axis xy; axis square; axis off;
title(sprintf('%0.2f', peak_rate), 'fontname', 'calibri', 'HorizontalAlignment', 'left')

subplot(3,4,10)
imagesc(rate_mat); axis xy; axis square; axis off;
caxis([0 lim]);
title(sprintf('%0.2f', lim), 'fontname', 'calibri', 'HorizontalAlignment', 'left')

%% 3rd column

subplot(3,4,3)
plot(pos_mean_x,pos_mean_y,'k'); axis image; hold on;
plot(spkx,spky,'.r');
axis off;
title('140 mins post-inj.', 'fontname', 'calibri')

subplot(3,4,7)
imagesc(rate_mat); axis xy; axis square; axis off;
title(sprintf('%0.2f', peak_rate), 'fontname', 'calibri', 'HorizontalAlignment', 'left')

subplot(3,4,11)
imagesc(rate_mat); axis xy; axis square; axis off;
caxis([0 lim]);
title(sprintf('%0.2f', lim), 'fontname', 'calibri', 'HorizontalAlignment', 'left')


%% 4th column

subplot(3,4,4)
plot(pos_mean_x,pos_mean_y,'k'); axis image; hold on;
plot(spkx,spky,'.r');
axis off;
title('recovery from CNO', 'fontname', 'calibri')

subplot(3,4,8)
imagesc(rate_mat); axis xy; axis square; axis off;
title(sprintf('%0.2f', peak_rate), 'fontname', 'calibri', 'HorizontalAlignment', 'left')

subplot(3,4,12)
imagesc(rate_mat); axis xy; axis square; axis off;
caxis([0 lim]);
title(sprintf('%0.2f', lim), 'fontname', 'calibri', 'HorizontalAlignment', 'left')