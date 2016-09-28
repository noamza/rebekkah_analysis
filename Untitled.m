function PlotShuffleHistogram(n,m,o, shuffle_data,real_value, y_lim,y_label,x_label, letter)

subplot(n,m,o)
hist(shuffle_data); hold on; box off;
y= 0:y_lim;
x= ones(1,y_lim+1)* real_value; 
plot(x,y,'-', 'LineWidth', 2);
ylabel(y_label)
xlabel(x_label)
text(-0.25,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
ylim([0 y_lim])
