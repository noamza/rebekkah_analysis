function fano=InverseFourierVariabilityAnalysis(inverse_mat)

autocorr= Cross_Correlation(inverse_mat,inverse_mat);
auto_inds= FindAutoMaxInds(autocorr);
PF_radius=findPlaceFieldRadius(autocorr, auto_inds);

fano= findFano(inverse_mat,PF_radius);



