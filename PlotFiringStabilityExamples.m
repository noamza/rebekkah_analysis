function PlotFiringStabilityExamples

% plot within session results
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF3')
load('info and results of half sessions G3MD25PF3.mat')
%best_inds=find(all_corrs>0.7 & all_corrs2>0.7 & twoD_zone_corrs>0.7);

figure;

best_inds=[64 213 192];

for ex_num=1:2
    add=(ex_num-1)*3;
    ImagePlot(ex_num,rate_mats_b,rate_mats_e,peak_rates_b,peak_rates_e,best_inds,add,'A')
end
clearvars

%plot same session results
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info')
load('corr coef results of same arenas updated.mat')

best_inds=[9 5];

for ex_num=1:2
    add=(ex_num+2-1)*3;
    ImagePlot(ex_num,all_rm_1,all_rm_2,all_rates_1,all_rates_2,best_inds,add, 'B')
end
clearvars

%plot rescaling session results
load('corr coef results of rescaled arenas.mat')

best_inds=find(all_corrs>0.9);
%best_inds=best_inds(3:end);

for ex_num=1:2
    add=(ex_num+4-1)*3;
    ImagePlot(ex_num,all_rm_1,all_rm_2,all_rates_1,all_rates_2,best_inds,add, 'C')
end


set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 30.45], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 30.45])


function ImagePlot(ex_num,rate_mats_b,rate_mats_e,peak_rates_b,peak_rates_e,best_inds,add,letter)

n=6;
m=3;

ind=best_inds(ex_num);

subplot(n,m,1+add)
imagesc(rate_mats_b{ind});
axis off;
axis square;

if add+1==1 || add+1==7 || add+1==13
text(-0.3,1.1,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
end

subplot(n,m,2+add)
imagesc(rate_mats_e{ind})
axis off;

if add+2==14 || add+2==17
   axis image
else
axis square;
end

if ex_num==1
    subplot(n,m,1+add)
    title('First half','HorizontalAlignment', 'right')
    
    subplot(n,m,2+add)
    title('Second half','HorizontalAlignment', 'right')
end

rates_1=peak_rates_b{ind};
[rates_1,rank]=sort(rates_1);
rates_2=peak_rates_e{ind};
rates_2=rates_2(rank);

max_lim=max([rates_1(end) rates_2(end)]);

subplot(n,m,3+add)
plot(rates_1,'o-','Color','black','MarkerFaceColor',[0.15 0.23 0.37]);hold on;
plot(rates_2,'o-','Color','black','MarkerFaceColor',[0.73 0.83 0.96]);
%xlabel('1st half')
ylabel('firing rate [Hz]')
ylim([0 max_lim])
xlim([1 length(rates_1)])
box off;


disp('')