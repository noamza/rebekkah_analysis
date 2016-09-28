function GaussianMatsAndPFRadii

load('3sets G3MD25 data and results.mat')

gaussian_mats=cell(1,length(fanos));
for i =1:length(fanos)
    
    pos_x=pos_x_all{i};
    pos_y=pos_y_all{i};
    pos_t=pos_t_all{i};
    spk_x=spk_x_all{i};
    spk_y=spk_y_all{i};
    spk_t=spk_t_all{i};
    
    max_inds=max_indices{i};
    rate_mat=rate_mats_all{i}; 
    peak_rates=peak_rates_all{i};
    
    gaussian_mats{i}= createGaussianMat(rate_mat, PF_radii(i), max_inds, mean(peak_rates));
 
end

save('G3MD25PF7 gaussian mats and PF rads', 'gaussian_mats','PF_radii')