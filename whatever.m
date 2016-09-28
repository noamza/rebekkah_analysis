function whatever
n=1;m=2;

load('simulated results 367cells 100x per cell 3.mat')
load('3sets G3MD25PF3 data and results.mat')

all_fano_factors=mean(all_fano_factors);

hold on;
subplot(n,m,2)
hold on;
 PlotShuffleHistogram(all_fano_factors, 'Mean Fano factor','Simulations',...
     'B',mean(fanos),3000)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 11.6], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 11.6])



function PlotHistogram(data, x_label, y_label,letter)
hist(data);
box off;
xlabel(x_label)
ylabel(y_label)
text(-0.25,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)

function PlotShuffleHistogram(data, x_label, y_label,letter,real_value,y_lim)
hist(data); hold on;
y=0:y_lim;
x=ones(1,y_lim+1)*real_value;
plot(x,y,'-')
xlabel(x_label)
ylabel(y_label)
ylim([0 y_lim])
box off;
text(-0.25,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)