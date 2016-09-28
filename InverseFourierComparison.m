function InverseFourierComparison 

cd('C:\Users\Dori\Desktop\New folder\Dropbox\Ismakov et al\images') 
load('reconstructed inverse fourier.mat')

for i=1:length(recon_rms)

    inverse_mat=recon_rms{i};
    inverse_orig_mat=recon_orig_rms{i};
    
    fanos(i)=InverseFourierVariabilityAnalysis(inverse_mat);
    fanos_orig(i)=InverseFourierVariabilityAnalysis(inverse_orig_mat);

end    

a=fanos;

save('fanos of reconstructed fourtrans', 'fanos','fanos_orig')

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat', 'fanos')

c=fanos;

means=[mean(c),  mean(a)];
stderrs=[std(c)/sqrt(length(c)), std(a)/sqrt(length(a))];

figure
subplot(2,2,3)
errorbar([1, 1.5],means,stderrs,'o'); hold on;
errorbar(mean(c),std(c)/sqrt(length(c)),'o')

subplot(2,2,[1 2])
OverlayingHists(a,c)



disp('')
