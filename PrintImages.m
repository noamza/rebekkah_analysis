function PrintImages

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
load('3sets G3MD25 data and results.mat')

fourier_mats=cell(1,length(zone_mats_all));
recon_rms=cell(1,length(zone_mats_all));
recon_orig_rms=cell(1,length(zone_mats_all));
fourier_grids=cell(1,length(zone_mats_all));
fanos2=nan(1,length(zone_mats_all));

for i=1:length(zone_mats_all)
    
    fig=figure;
    n=2;m=5;
    
    % plot spike plot
    subplot(n,m,1)
    plot(pos_x_all{i},pos_y_all{i},'k');hold on;
    plot(spk_x_all{i},spk_y_all{i},'.r');
    axis equal;axis off;
    axis ij;
    title(sprintf('cell %d', i));
    
    %plot rate mat
    subplot(n,m,2)
    imagesc(rate_mats_all{i});
    axis equal; axis off;
    title(sprintf('%0.1f Hz', max(peak_rates_sm_all{i})), 'HorizontalAlignment', 'left');
    
    %plot zone mat
    subplot(n,m,3)
    imagesc(zone_mats_all{i});
    title('Zone Matrix')
    axis equal; axis off;
    
    %plot sorted rates
    subplot(n,m,6)
    plot(1:length(peak_rates_sm_all{i}),sort(peak_rates_sm_all{i}), 'ko-');
    title(sprintf('Fano factor=%0.1f', fanos(i)));
    axis square; box off;
    
    %plot autocorr
    subplot(n,m,7)
    imagesc(autocorrs_all{i});
    axis equal; axis off; axis ij;
    title('Autocorrelation')
    
    % plot fourier mat
    subplot(n,m,8)
    rate_mat=rate_mats_ns_all{i};
    mean_firing=nanmean(rate_mat(:));
    rate_mat=rate_mat-mean_firing;
    
    rate_mat(isnan(rate_mat))=0;
    
    % change rate map to size of 256 by 256 with zero padding
    size_rm= size(rate_mat);
    add_x=round((256-size_rm(1))/2);
    add_y=round((256-size_rm(2))/2);
    
    rate_mat_new=zeros(256,256);
    rate_mat_new(add_x:add_x+size_rm(1)-1,add_y:add_y+size_rm(2)-1)=rate_mat;
    
    rate_mat=rate_mat_new;
    
    % create fourier matrix
    fourier_mat=fft2(rate_mat);
    
    f_m_orig=fourier_mat;
    
    fourier_mat(1)=0;
    [x,y] = size(fourier_mat);
    fourier_mat = fftshift(fourier_mat);
    x1 = round(5*x/20); x2 = round(15*x/20); y1 = round(5*y/20); y2 = round(15*y/20);
    
    fm_show=fourier_mat(x1:x2,y1:y2);
    imagesc(abs(fm_show));
    axis equal;axis off; axis ij;
    
    % create simplified fourier mat
    subplot(n,m,9)
    
    fourier_mat_new=CreateFourierTransformMat(fourier_mat,x1,x2,y1,y2);
    
    %max_f= max(abs(fourier_mat(:)));
    %fourier_mat_new=fourier_mat;
    %fourier_mat_new(abs(fourier_mat)<max_f*0.7)=0;
    fm_show_new=fourier_mat_new(x1:x2,y1:y2);
    imagesc(abs(fm_show_new));
    axis image; axis off; axis ij;
    
    fourier_mats{i}=fm_show;
    fourier_grids{i}=fm_show_new;
    
    subplot(n,m,4)
    inverse_fft=ifft2(f_m_orig);
    inverse_fft(1)=f_m_orig(1);
    
    inverse_fft=inverse_fft(add_x:add_x+size_rm(1)-1,add_y:add_y+size_rm(2)-1);
     parms.sigma=1;
     inverse_fft=SmoothRateMat(inverse_fft,parms);
    
    
    inverse_fft= StretchImage(inverse_fft, size(inverse_fft), size(inverse_fft)*2);
    
    imagesc((inverse_fft));
    axis equal; axis off; axis ij;
    title('Reconstucted orig')
    
    recon_orig_rms{i}=inverse_fft;
    
    subplot(n,m,5)
    
    inverse_fft=ifftshift(fourier_mat_new);
    inverse_fft(1)=f_m_orig(1);
    inverse_fft=ifft2(inverse_fft);
    %inverse_fft=inverse_fft(add_x:add_x+size_rm(1)-1,add_y:add_y+size_rm(2)-1);
    inverse_fft=inverse_fft(add_x:add_x+size_rm(1)-1,add_y:add_y+size_rm(2)-1);
    parms.sigma=1;
    inverse_fft=SmoothRateMat(inverse_fft,parms);
    
    
    inverse_fft= StretchImage(inverse_fft, size(inverse_fft), size(inverse_fft)*2);
    
    imagesc((inverse_fft));
    axis equal; axis off; axis ij;
    title('Reconstucted edited')
    
    %check example, make sure (1) is largest value. check if not
    
    recon_rms{i}=inverse_fft;
   
    fano=InverseFourierVariabilityAnalysis(inverse_fft);
    max_inds=FindMaxIndsRateMap(inverse_fft);
    peak_rates=findPeakRates(max_inds, inverse_fft);
    fanos2(i)=fano;
    
     %plot sorted rates
    subplot(n,m,10)
    plot(1:length(peak_rates),sort(peak_rates), 'ko-');
    title(sprintf('Fano factor=%0.1f', fano));
    axis square; box off;
    
    %save images
    cd('C:\Users\Dori\Desktop\New folder\Dropbox\Ismakov et al\images');
    saveas(fig,sprintf('%d.jpg',i));
    
    i
    
    close all
end

save('reconstructed inverse fourier', 'fanos2','recon_rms', 'recon_orig_rms', 'fourier_grids', 'fourier_mats')

