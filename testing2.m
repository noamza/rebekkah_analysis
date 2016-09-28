function testing 

load('remapping vs noremapping corrs.mat')

    a=no_remap_corrs;
    b=remap_corrs;
    
    n=2;m=3;
    figure;
    subplot(n,m,1)
    PlotMeanErrorbars(a,b);

bin_values=0:0.2:1.1;
    for h=1:length(bin_values)-1
        x=xcorr_scores;
    bins=find(x>=bin_values(h) &x<bin_values(h+1));
    binned_data(h)=mean(all_corrcoefs(bins));
    stderrs(h)=FindStdErr(all_corrcoefs(bins));
    mean_x(h)=mean(xcorr_scores(bins));
    end
    
   subplot(n,m,[2 3])
   errorbar(mean_x,binned_data,stderrs);
   
   subplot(n,m,[4 5])
  OverlayingHists(a,b)
    

    
    function PlotMeanErrorbars(a,b)
    
    all_means= [mean(a) mean(b)];
    stderr1=FindStdErr(a);
    stderr2=FindStdErr(b);
    %stderr3=FindStdErr(c);
    all_stderrs= [stderr1 stderr2];
    
    errorbar(all_means, all_stderrs, 'o')
    
        function stderr= FindStdErr(x)
            
            stderr= nanstd(x)/ sqrt((sum(~isnan(x))));
            
            