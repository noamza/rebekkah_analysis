function GaussianMatFromRateMat 

load('3sets G3MD15PF7 data and results.mat', ...
    'max_indices', 'rate_mats_all', 'peak_rates_all',...
    'PF_radii')

gaussian_mats=cell(1,length(max_indices));
for i =1:length(max_indices)
    
    max_inds=max_indices{i};
    rate_mat=rate_mats_all{i}; 
    peak_rates=peak_rates_all{i};
    
    gaussian_mats{i}= createGaussianMat(rate_mat, PF_radii(i), max_inds, mean(peak_rates));
 
end

save('G3MD15PF7 gaussian mats', 'gaussian_mats')