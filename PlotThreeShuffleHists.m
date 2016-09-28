function PlotThreeShuffleHists(hist_1,hist_2,hist_3, ...
    real_value, real_value2,real_value3,...
    y_lim, y_label, x_label, y_label2, x_label2, y_label3, x_label3) 

figure;
n=3;
m=2

subplot(n,m,1)
hist(hist_1); hold on; box off;
y= 0:y_lim;
x= ones(1,y_lim+1)* real_value; 
plot(x,y,'-', 'LineWidth', 2);
ylabel(y_label)
xlabel(x_label)
text(-0.25,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
ylim([0 y_lim])

subplot(n,m,3)
hist(hist_2); hold on; box off;
y= 0:y_lim;
x= ones(1,y_lim+1)* real_value2; 
plot(x,y,'-', 'LineWidth', 2);
ylabel(y_label2)
xlabel(x_label2)
text(-0.25,1.05,'B','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
ylim([0 y_lim])

subplot(n,m,5)
hist(hist_3); hold on; box off;
y= 0:y_lim;
x= ones(1,y_lim+1)* real_value3;
plot(x,y,'-', 'LineWidth', 2);
ylabel(y_label3);
xlabel(x_label3); 
text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
ylim([0 y_lim])





