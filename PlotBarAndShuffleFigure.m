function PlotBarAndShuffleFigure(bar_graph,bar_x,shuffle_hist, ...
    real_value, ...
    y_lim, y_label, x_label, y_label2, x_label2) 

figure;
n=1;
m=2;

subplot(n,m,1)
bar(bar_x,bar_graph); hold on; box off;
ylabel(y_label)
xlabel(x_label)
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

subplot(n,m,2)
hist(shuffle_hist); hold on; box off;
y= 0:y_lim;
x= ones(1,y_lim+1)* real_value; 
plot(x,y,'-', 'LineWidth', 2);
ylabel(y_label2)
xlabel(x_label2)
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
ylim([0 y_lim])

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 5.8], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 5.8])
