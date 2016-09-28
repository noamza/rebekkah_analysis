%Plot variability plot

cd('N:\users\rebekkah\results and info of analysis')
load('norm variability analysis results.mat')

figure;

d=errorbar(order, means_by_bins, stderr, 'k.');
hold on;
box off;
xlabel('Normalized field rank', 'FontSize', 20)
ylabel('Normalized firing rate', 'FontSize', 20)
 

all_norm_orders= all_norm_orders';
all_sorted_rates= all_sorted_rates';
f= fit(all_norm_orders,all_sorted_rates,'poly1');
p1= 0.7748;
p2= 0.0674;
p1_1= 0.7423;
p2_1= 0.04754;
p1_2= 0.8074;
p2_2= 0.08726;

x=0:0.1:1.1;

y=(p1*x +p2);
y1= (p1_1*x +p2_1);
y2= (p1_2*x +p2_2);

h= plot(x,y, 'Color', [0.85 0.16 0], 'LineWidth', 2);
h1= plot(x,y1,'--', 'Color', [0.85 0.16 0],'LineWidth', 2)
h2= plot(x,y2,'--','Color', [0.85 0.16 0],'LineWidth', 2);
text(-0.2,1.05,'A','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 28)

legend([d h h1],{'mean firing rates','linear fit', '95% conf. intervals'}, 'FontSize', 20);

%set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 11.6, 8.4], 'PaperUnits', 'centimeters', 'PaperSize', [11.6, 8.4])

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 34.8, 23.2], 'PaperUnits', 'centimeters', 'PaperSize', [34.8, 23.2])