pre_inj_peak_rate= [7.28 6.36 2.24 1.93 6.95 4.95 3.63 2.69];
post_inj_peak_rate= [1.08 5.26 1.61 0.98 2.00 3.33 0.56 1.13];
recovery_peak_rate= [9.21 5.53 2.18 1.10 7.50 9.37 3.99 5.63];


pre_inj_mean_rate= [1.90 1.05 0.44 0.67 0.92 1.37 0.62 0.53];
post_inj_mean_rate= [0.08 0.70 0.40 0.41 0.16 0.58 0.12 0.27];
recovery_mean_rate=[1.36 0.89 0.47 0.40 0.50 1.50 0.61 1.80];

num= 1:8;


for num=1:8;
    
    y=[pre_inj_mean_rate(num) post_inj_mean_rate(num) recovery_mean_rate(num)] ;
    x=[ 1 2 3];
    
    plot(x,y, '-ok'); hold on;
    
end
    


