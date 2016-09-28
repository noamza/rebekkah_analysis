%use 25 4 threshs


context_means= [0.3773 0.2267 0.1556];
context_stderrs=[ 0.0305 0.0325 0.0305];  
category= {'SS', 'SC', 'D'};

n=2;m=2;

figure;
subplot(n,m,1)
errorbar(context_means, context_stderrs,'o', 'color', [0 0.75 0.75]);
box off;
ylabel('2D correlation');
xlim([0.7 3.3]); 
text(-0.25,1.05, 'A', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
change={'Diff Color'; 'Diff Scent'; 'Both Diff'}; 
set(gca,'XTickLabel',change,'XTick',1:3)


subplot(n,m,2)
errorbar(mean_ranks, stderrs, 'o', 'color', [0.85 0.16 0]);
box off;
ylabel('Distance rank');
xlim([0.7 3.3]); 
text(-0.25,1.05, 'B', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
change={'Same Context'; 'No Remap'; 'Remap'};
set(gca,'XTickLabel',change,'XTick',1:3)

subplot(n,m,3)
load('corrs pt 25 4 FIXED.mat')
errorbar(means, stderrs, 'o' , 'color', [0.85 0.16 0]);
box off;
ylabel('Corr coefficient');
xlim([0.7 3.3]); 
text(-0.25,1.05, 'C', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
set(gca,'XTickLabel',change,'XTick',1:3)

subplot(n,m,4)
load('binned corrs to mean dists UPDATED.mat')
errorbar(means,stderrs, 'ko-');
box off;
xlim([0.7 4.3]);
ylabel('Distance');
xlabel('2D correlation');
text(-0.25,1.05, 'D', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)




load('dists and corrs 15 25 2 all.mat')
load('binned data.mat')
n=2;m=2;
figure;

change= [{'S.C.'} {'N.R.'} {'R'}];

subplot(n,m,1)
errorbar(1, corr_means(1), corr_stderrs(1), 'o', 'color', [0.75 0.75 0]); 
hold on;
errorbar(2, corr_means(2), corr_stderrs(2), 'o','color', [0.49 0.49 0.49]);
errorbar(3, corr_means(3), corr_stderrs(3), 'o','color', [0.51 0.38 0.48]);
box off;
ylabel('Corr coefficient');
xlim([0.7 3.3]);

text(-0.25,1.05, 'A', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
set(gca,'XTickLabel',change,'XTick',1:3)
legend('Same Context', 'Non-Remapped', 'Remapped')
legend('boxoff')


subplot(n,m,2)
errorbar(lims(1:7), means,stderrs, 'ko-');
box off;
xlim([lims(1)-0.1  lims(7)+0.1]);
ylim([-0.25 1])
ylabel('Corr coeff');
xlabel('2D correlation');
text(-0.25,1.05, 'B', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
legend('Binned data')
legend('boxoff')

subplot(n,m,3)
hist(no_remap_corrs); hold on;
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 0.75 0.75],'EdgeColor','w','facealpha',0.75)
hold on;
hist(remap_corrs)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75, 'FaceColor', [0.49 0.49 0.49]);
legend('Non-Remapped', 'Remapped')
legend('boxoff')

ylim([0 65]);

text(-0.25,1.05, 'C', 'Units', 'Normalized', ...
    'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 11.6], ...
    'PaperUnits', 'centimeters', 'PaperSize', [17.4 11.6])