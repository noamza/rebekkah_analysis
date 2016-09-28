function PlotShuffleHistogram(n,m,o, shuffle_data,real_value, y_lim,y_label,x_label, letter)

subplot(n,m,o)
hist(shuffle_data, 10, [], [0.76 0.87 0.78]); hold on; box off;
y= 0:y_lim;
x= ones(1,y_lim+1)* real_value; 
plot(x,y,'-', 'LineWidth', 2, 'color', [0.5 0.5 0.5]);
ylabel(y_label)
xlabel(x_label)
text(-0.25,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
ylim([0 y_lim])
