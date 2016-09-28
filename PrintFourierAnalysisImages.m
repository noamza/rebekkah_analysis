function PrintFourierAnalysisImages

dbstop if error


cd('C:\Users\Dori\Desktop\New folder\Dropbox\Ismakov et al\images')
A=load('reconstructed inverse fourier.mat', 'fanos', 'fourier_grids',...
    'recon_rms', 'fourier_mats');

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
B=load('3sets G3MD25 data and results.mat', 'fanos', ...
    'peak_rates_sm_all','rate_mats_all', 'norm_max_index', 'norm_dist');

best_examples=[14,25,30]; 

figure;

n=6; m=4;
letters={{'A'},{'B'},{'C'}};
add=0;

step=8;

for h=1:3

    ind=best_examples(h); 
    
    letter=letters{h}; 
    rm= B.rate_mats_all{ind}; 
    pr= B.peak_rates_sm_all{ind}; 
    ft= abs(A.fourier_grids{ind});
    rrm= A.recon_rms{ind};
    f= abs(A.fourier_mats{ind}); 
    
    % plot rate map
    subplot(n,m,[1]+(add*step));
    imagesc(rm);
    axis image; axis off; 
    
    text(-0.45,1.05,letter,'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    
    % plot sorted rates
    subplot(n,m,[5 6]+ (add*step));
    plot(1:length(pr),sort(pr), 'ko-', 'MarkerFaceColor',[0.76,0.87,0.78]);
    ylabel('Firing rate [Hz]') 
    title(sprintf('Fano factor=%0.1f',B.fanos(ind))); 
    box off; 
    xlim([1 length(pr)])
    ylim([0 max(pr)])
    % plot fourier mat
    subplot(n,m,[2]+(add*step));
    imagesc(f);
    axis image; axis off; 
    % plot fourier grid components
    subplot(n,m,[3]+(add*step));
    imagesc(ft);
    axis image; axis off; 
    % plot reconstructed rate map
    subplot(n,m,[4]+(add*step));
    imagesc(rrm);
    axis image; axis off; 
    % plot reconstructed sorted rates
    max_inds=FindMaxIndsRateMap(rrm);
    rpr= findPeakRates(max_inds, rrm);
    
    subplot(n,m,[7 8]+ (add*step));
    plot(1:length(rpr),sort(rpr), 'ko-');
    ylabel('Firing rate [Hz]') 
    title(sprintf('Fano factor=%0.1f',A.fanos(ind))); 
    box off; 
    xlim([1 length(rpr)])
    ylim([0 max(pr)])
    
    add=add+1;
end 

figure; 
subplot(4,6,[1 2 9 10])
a=A.fanos;
b=B.fanos;
% 
 means=[mean(b) mean(a)] ;
 stderrs= [std(b)/sqrt(length(b))  std(a)/sqrt(length(a))];
% 
 errorbar(means, stderrs,'ko'); hold on; 
 errorbar(mean(b), stderrs(1),'ko', 'MarkerFaceColor',[0.76,0.87,0.78]); hold on; 
 box off;
% 
  text(-0.45,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 23.2], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 23.3])


