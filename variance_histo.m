function [ output_args ] = Untitled6( input_args )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here



% get all files in folder
% create file_list
% create loop 

parms.dir_load_data = 'C:\Users\Dori\Desktop\Rebekkah data\save results with where'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};


% load specific file.
% from each fileget S.variance. store in vector variance.

% figure;
% 
% n=1;
% m=3; 
% 

for i=1:length(file_names)
    file_name = file_names{i};
    load(file_name);

      variances(i)= S.variance;    
%     gridness(i)= S.gridness2;
      rayleigh_HD(i) = S.HD.rayleigh.score;
%     rayleigh_MD(i)= S.MD.rayleigh.score;


peak_means(i) = mean(S.peak_rates);

% %      
% %   
% 
% 
% 
% 
% radii (i)= S.PF_radius;
% 

% 
% mean_norm_means(i) = mean(S.norm_means);
% 
% variance_norm_means(i)= var(S.norm_means);
% 
% mean_means(i) = mean(S.sorted_means);   
%     
% % min_mean_max(1)= min(S.norm_means);
% % min_mean_max(2)= mean(S.norm_means);
% % min_mean_max (3) = max(S.norm_means);
% 
% % plot(min_mean_max, 'o-', 'color', 'g');
% % hold on;
% % drawnow;
% 
% % min_mean_max(1)= min(S.sorted_means);
% % min_mean_max(2)= mean(S.sorted_means);
% % min_min_max (3) = max(S.sorted_means);
% 
% %     
% %     
% %     
%      if rayleigh_HD (i) > 0.3
%         plot(min_mean_max, 'o-','color','r');
%        hold on;
%      elseif rayleigh_HD (i) < 0.3
%           plot(min_mean_max, 'o-', 'color','b');
%           hold on;
%      else
%         disp ('');
%      end
% 
%     drawnow;
% 
%     
% %     
% %     S.mean = S.mean'
% %     avg_mean = mean(S.mean)
% %     means(i) = avg_mean;
% %     
% %     clear S;
% 
    
locations(i) = S.where;
% 

end

%
% disp('');
% 
% std = sqrt(variances2); 
% 
% 
for i= 1:220;
% 
var_over_mean(i) = variances(i)/peak_means(i);
% 
% std_mean (i) = std(i)/mean_means(i);
% 
end
% 
% 
log_var_mean = log10(var_over_mean);



%hist(log_var_mean)


i= 1:220;
% 
% no_HD_inds = find(rayleigh_HD(i) > 0.3 ) ;
% no_HD_norm_values(i) = mean_norm_means_list(i);
% no_HD_norm_values(no_HD_inds)= [];
% no_HD_norm_values(no_HD_norm_values==0) =[];
% 
% HD_inds = find(rayleigh_HD(i) < 0.3 ) ;
% HD_norm_values(i) = mean_norm_means_list(i);
% HD_norm_values(HD_inds)= [];
% HD_norm_values(HD_norm_values==0) =[];


no_HD_inds = find(rayleigh_HD(i) > 0.3 ) ;
HD_inds = find(rayleigh_HD(i) < 0.3 ) ;
low_v_o_m_inds = find(var_over_mean(i) > 2);
high_v_o_m_inds = find(var_over_mean(i) < 2);
no_HD_low_var_inds = union(no_HD_inds, low_v_o_m_inds);
no_HD_high_var_inds = union (no_HD_inds, high_v_o_m_inds);
HD_high_var_inds = union (HD_inds, high_v_o_m_inds);
HD_low_var_inds= union (HD_inds, low_v_o_m_inds);


var_over_mean = log_var_mean;

%pure grid cells of all variability

no_HD_all_var(i) = var_over_mean(i);
no_HD_all_var(no_HD_inds)= [];
no_HD_all_var(no_HD_all_var==0) =[];

%conjunctive gric cells of all variability
HD_all_var(i) = var_over_mean(i);
HD_all_var(HD_inds) = [];
HD_all_var(HD_all_var==0) =[];

%low variability cell of all directionality scores
low_var_all_grid = var_over_mean(i);
low_var_all_grid(low_v_o_m_inds) = [];
low_var_all_grid(low_var_all_grid==0) =[];

%high variability cell of all directionality scores
high_var_all_grid = var_over_mean(i);
high_var_all_grid(high_v_o_m_inds) = [];
high_var_all_grid(high_var_all_grid==0) =[];

%pure grid cells of high variability 
no_HD_high_var = var_over_mean(i);
no_HD_high_var (no_HD_high_var_inds) = [];
no_HD_high_var(no_HD_high_var==0)=[];

%pure grid cells of low variability
no_HD_low_var = var_over_mean(i);
no_HD_low_var (no_HD_low_var_inds) = [];
no_HD_low_var(no_HD_low_var==0)=[];

%conjunctive grid cells of high variability
HD_high_var = var_over_mean(i);
HD_high_var (HD_high_var_inds) = [];
HD_high_var(HD_high_var==0)=[];

%conjunctive grid cells of low variability
HD_low_var = var_over_mean(i);
HD_low_var (HD_low_var_inds) = [];
HD_low_var(HD_low_var==0)=[];

% 
% %same with std_mean
% 
% % std_no_HD_var(i) = std_mean(i);
% % std_no_HD_var(no_HD_var_inds)= [];
% % std_no_HD_var(no_HD_var==0) =[];
% % 
% % std_HD_var(i) = std_mean(i);
% % std_HD_var(HD_var_inds)= [];
% % std_HD_var(HD_var==0) =[];
% 
% 
figure;

n=3;
m=3;

subplot (n,m,1)
hist(no_HD_low_var(:));
title('pure grid cells of low variability')

subplot (n,m,2)
hist(no_HD_high_var(:));
title('pure grid cells of high variability')

subplot (n,m,3)
hist(no_HD_all_var(:));
title('pure grid cells')

subplot (n,m,4)
hist(HD_low_var(:));
title('conjunctive grid cells of low variability')

subplot (n,m,5)
hist(HD_high_var(:));
title('conjunctive grid cells of high variability')

subplot (n,m,6)
hist(HD_all_var(:));
title('conjunctive grid cells')

subplot (n,m,7)
hist(low_var_all_grid(:));
title('grid cells of low variability')

subplot (n,m,8)
hist(high_var_all_grid(:));
title('grid cells of high variability')

subplot (n,m,9)
hist(var_over_mean(:));
title('total grid cell variability')


% subplot(n,m,1)
% xvalues = 0:0.02:1;
% hist(no_HD_norm_values(:), xvalues);
% axis ([0 1 0 35]);
% title('means of normalized mean firing rate values of grid cells with low directionality');
%    
% subplot(n,m,2)
% xvalues = 0:0.02:1;
% hist(HD_norm_values(:), xvalues);
% axis ([0 1 0 35]);
% title('means of normalized mean firing rate values of grid cells with high directionality');
% 
% subplot(n,m,3)
% xvalues = 0:0.02:1;
% hist(mean_norm_means_list(:), xvalues);
% axis ([0 1 0 35]);
% title('means of normalized mean firing rate values of grid cells'); 
% 
% % 
% % figure;
% % n=4
% % m=1
% % 
% % subplot(n,m,1) 
% % %xvalues= 0:25;
% % hist(no_HD_var(:),-2:0.15:2);
% % hold on;
% %  
% % %xvalues= 0:25; 
% % hist(HD_var(:), -2:0.15:2);
% % 
% % h = findobj(gca, 'Type','patch');
% % set(h(1), 'FaceColor','r', 'EdgeColor','w')
% % set(h(2), 'FaceColor','b', 'EdgeColor','w')
% % title('variance over mean in grid cells'); 
% 
% % subplot(n,m,2) 
% % %xvalues= 0:0.025:2;
% % hist(std_no_HD_var(:), xvalues);
% % hold on;
% %  
% % xvalues= 0:0.025:2; 
% % hist(std_HD_var(:), xvalues);
% % 
% % h = findobj(gca, 'Type','patch');
% % set(h(1), 'FaceColor','r', 'EdgeColor','w')
% % set(h(2), 'FaceColor','b', 'EdgeColor','w')
% % title('standard deviation over mean in grid cells');
% % 
% 
% subplot(n,m,3) 
% xvalues= 0:0.05:2;
% hist(std_no_HD_var(:), xvalues);
% hold on;
%  
% xvalues= 0:0.05:2; 
% hist(std_HD_var(:), xvalues);
% 
% h = findobj(gca, 'Type','patch');
% set(h(1), 'FaceColor','r', 'EdgeColor','w')
% set(h(2), 'FaceColor','b', 'EdgeColor','w')
% title('standard deviation over mean in grid cells');
% 
% 
% subplot(n,m,4) 
% xvalues= 0:0.1:2;
% hist(std_no_HD_var(:), xvalues);
% hold on;
%  
% xvalues= 0:0.1:2; 
% hist(std_HD_var(:), xvalues);
% 
% h = findobj(gca, 'Type','patch');
% set(h(1), 'FaceColor','r', 'EdgeColor','w')
% set(h(2), 'FaceColor','b', 'EdgeColor','w')
% title('standard deviation over mean in grid cells');

%subplot(n,m,2)
% xvalues= 0:5:350; 
% hist(var_over_mean(:), xvalues);
% %axis ([0 1 0 35]);
% title('variablity in grid cells'); 





% low_grid_var(i) = variances(i)
% low_gridness(i) = gridness(i)
% 
% low_grid_var(low_inds)= [];
% low_grid_var(low_grid_var==0) = [];
% low_gridness(low_inds) = [];
% low_gridness(low_gridness==0)=[];
% 
% med_inds = find(gridness(i)<0.2 | gridness(i)>0.6 | isnan(gridness(i))) ;
% med_grid_var(i) = variances(i)
% med_grid_var(med_inds) = [];
% med_grid_var(med_grid_var==0) = [];
% 
% hi_med_inds = find(gridness(i)<0.6 | gridness(i)>0.9 | isnan(gridness(i))) ;
% hi_med_grid_var(i) = variances(i)
% hi_med_grid_var(hi_med_inds) = [];
% hi_med_grid_var(hi_med_grid_var==0) = [];
% 


%looking at variances seperated into 3 group of gridness levels

% i=1:743
% 
% low_inds = find(gridness(i)> 0.2 | gridness(i) < 0 | isnan(gridness(i)) )
% 
% low_grid_var(i) = variances(i)
% low_gridness(i) = gridness(i)
% 
% low_grid_var(low_inds)= [];
% low_grid_var(low_grid_var==0) = [];
% low_gridness(low_inds) = [];
% low_gridness(low_gridness==0)=[];
% 
% med_inds = find(gridness(i)<0.2 | gridness(i)>0.6 | isnan(gridness(i))) ;
% med_grid_var(i) = variances(i)
% med_grid_var(med_inds) = [];
% med_grid_var(med_grid_var==0) = [];
% 
% hi_med_inds = find(gridness(i)<0.6 | gridness(i)>0.9 | isnan(gridness(i))) ;
% hi_med_grid_var(i) = variances(i)
% hi_med_grid_var(hi_med_inds) = [];
% hi_med_grid_var(hi_med_grid_var==0) = [];
% 
% high_inds = find(gridness(i)<0.9 | isnan(gridness(i)));
% high_grid_var(i) = variances(i)
% high_grid_var(high_inds) = [];
% high_grid_var(high_grid_var==0) = [];
% 
% 
% low_inds = find(rayleigh_HD(i)>0.1 | isnan(rayleigh_HD(i)));
% low_HD_var(i) = variances(i)
% low_HD_var(low_inds) = [];
% low_HD_var(low_HD_var==0) = [];
% 
% med_inds = find(rayleigh_HD(i)<0.1 | rayleigh_HD(i)>0.2 | isnan(rayleigh_HD(i)));
% med_HD_var(i) = variances(i)
% med_HD_var(med_inds) = [];
% med_HD_var(med_HD_var==0) = [];
% 
% hi_med_inds = find(gridness(i)<0.2 | gridness(i)>0.5 | isnan(gridness(i))) ;
% hi_med_HD_var(i) = variances(i)
% hi_med_HD_var(hi_med_inds) = [];
% hi_med_HD_var(hi_med_HD_var==0) = [];
% 
% high_inds = find(rayleigh_HD(i)<0.5 | isnan(rayleigh_HD(i)));
% high_HD_var(i) = variances(i)
% high_HD_var(high_inds) = [];
% high_HD_var(high_HD_var==0) = [];
% 
% 
% low_inds = find(rayleigh_MD(i)>0.05 | isnan(rayleigh_MD(i)));
% low_MD_var(i) = variances(i)
% low_MD_var(low_inds) = [];
% low_MD_var(low_MD_var==0) = [];
% 
% med_inds = find(rayleigh_MD(i)<0.05 | rayleigh_MD(i)>0.1 | isnan(rayleigh_MD(i)));
% med_MD_var(i) = variances(i)
% med_MD_var(med_inds) = [];
% med_MD_var(med_MD_var==0) = [];
% 
% hi_med_inds = find(gridness(i)<0.1 | gridness(i)>0.4 | isnan(gridness(i))) ;
% hi_med_MD_var(i) = variances(i)
% hi_med_MD_var(hi_med_inds) = [];
% hi_med_MD_var(hi_med_MD_var==0) = [];
% 
% high_inds = find(rayleigh_MD(i)<0.4 | isnan(rayleigh_MD(i)));
% high_MD_var(i) = variances(i)
% high_MD_var(high_inds) = [];
% high_MD_var(high_MD_var==0) = [];
% 
% 
% n= 3;
% m= 4; 
% 
% max_lim = max(low_grid_var) + 10
% 
% figure;
%     subplot(n,m,1)
%     xvalues = 0:2.5:250;
%     hist(low_grid_var(:), xvalues);
%     xlim([0 max_lim])
%     title('variablilty of low gridness');
% 
%    max_lim = max(med_grid_var) + 10
%     
%     subplot(n,m,2)
%     xvalues = 0:2.5:250;
%     hist (med_grid_var(:), xvalues);
%    xlim([0 max_lim])
%     title ('variability of average gridness')
%    
%     
%       max_lim = max(hi_med_grid_var) + 10
%     
%     subplot(n,m,3)
%     xvalues = 0:2.5:250;
%     hist (hi_med_grid_var(:), xvalues);
%    xlim([0 max_lim])
%     title ('variability of med to high gridness')
%    
%     
%     
%     max_lim = max(high_grid_var) + 10
%     
%     subplot(n,m,4)
%     xvalues = 0:2.5:250;
%     hist(high_grid_var(:), xvalues);
%     xlim([0 max_lim])
%     title('variability of high gridness')
%       
%     max_lim = max(low_HD_var) + 10
%     
%     subplot(n,m,5)
%     xvalues = 0:2.5:250;
%     hist(low_HD_var(:), xvalues);
%    xlim([0 max_lim])
%     title('variablilty of low HD');
% 
%     max_lim = max(med_HD_var) + 10
%     
%     subplot(n,m,6)
%     xvalues = 0:2.5:250;
%     hist (med_HD_var(:), xvalues);
%     xlim([0 max_lim])
%     title ('variability of average HD')
%   
%     
%      max_lim = max(hi_med_HD_var) + 10
%     
%     subplot(n,m,7)
%     xvalues = 0:2.5:250;
%     hist (hi_med_HD_var(:), xvalues);
%    xlim([0 max_lim])
%     title ('variability of med to high HD')
%     
%     max_lim = max(high_HD_var) + 10
%     
%     subplot(n,m,8)
%     xvalues = 0:2.5:250;
%     hist(high_HD_var(:), xvalues);
%     xlim([0 max_lim])
%     title('variability of high HD')
%     
%     max_lim = max(low_MD_var) + 10
%     
%     subplot(n,m,9)
%     xvalues = 0:2.5:250;
%     hist(low_MD_var(:), xvalues);
%     xlim([0 max_lim])
%     title('variablilty of low MD');
% 
%     max_lim = max(med_MD_var) + 10
%     
%     subplot(n,m,10)
%     xvalues = 0:2.5:250;
%     hist (med_MD_var(:), xvalues);
%    xlim([0 max_lim])
%     title ('variability of average MD')
%    
%     max_lim = max(hi_med_HD_var) + 10
%     
%     subplot(n,m,11)
%     xvalues = 0:2.5:250;
%     hist (hi_med_MD_var(:), xvalues);
%    xlim([0 max_lim])
%     title ('variability of med to high MD')
%     
%     max_lim = max(high_MD_var) + 10
%     
%     subplot(n,m,12)
%     xvalues = 0:2.5:250;
%     hist(high_MD_var(:), xvalues);
%     xlim([0 max_lim])
%     title('variability of high MD')
% %    
    
    
% high_variance = find(variances > 50)
% 
% 
% histology_high_var = histology_layer(high_variance)  



% histogram of variability in all cells

% figure;
% xvalues = 0:2.5:100; 
% hist(variances(:), xvalues);
% xlim([0.1 100])
% title('Variability of firing rates between place fields of individual grid cells'); 
% % 
% % scatterplot of gridness vs. variablity
% 
% figure;
% scatter(gridness, rayleigh_HD)
% title('HD vs gridness number')
% lsline
% %scatterplot of HD vs variablity
% 
% figure;
% scatter(variances, rayleigh_HD)
% xlim([0.25 100])
% title ('variability vs HD rayleigh score')
% lsline
% lsqcurvefit
% 
% %scatterplot of MD vs variablity
% 
% figure;
% scatter(variances, rayleigh_MD)
% xlim([0.25 100])
% title ('variability vs MD rayleigh score')
% lsline
% 
% h= lsline 
% 
% lsqcurvefit(h, x0, variances, rayleigh_MD) 


% figure;
% scatter(variances, means)
% xlim([1 100])
% title('variability vs average mean of PF firing rates')
% lsline


disp('finished')

end

