
function AppendingResults

% A=load('same arenas corr coef results of rescaling data.mat');
% B=load('same arenas corr coef results.mat');
% 
% all_rates_1=AppendTwoVectors(A.all_rates_1,B.all_rates_1);
% all_rates_2=AppendTwoVectors(A.all_rates_2,B.all_rates_2);
% all_zm_1=AppendTwoVectors(A.all_zm_1,B.all_zm_1);
% all_zm_2=AppendTwoVectors(A.all_zm_2,B.all_zm_2);
% all_rm_1=AppendTwoVectors(A.all_rm_1,B.all_rm_1);
% all_rm_2=AppendTwoVectors(A.all_rm_2,B.all_rm_2);
% all_corrs=AppendTwoVectors(A.all_corrs,B.all_corrs);
% all_corrs2=AppendTwoVectors(A.all_corrs2,B.all_corrs2);
% all_max_inds_1=AppendTwoVectors(A.all_max_inds_1,B.all_max_inds_1);
% all_max_inds_2=AppendTwoVectors(A.all_max_inds_2,B.all_max_inds_2);
% orig_rates_1=AppendTwoVectors(A.orig_rates_1,B.orig_rates_1);
% orig_rates_2=AppendTwoVectors(A.orig_rates_2,B.orig_rates_2);
% 
% pos_x_all1=AppendTwoVectors(A.all_pos_x1,B.all_pos_x1);
% pos_x_all2=AppendTwoVectors(A.all_pos_x2,B.all_pos_x2);
% pos_y_all1=AppendTwoVectors(A.all_pos_y1,B.all_pos_y1);
% pos_y_all2=AppendTwoVectors(A.all_pos_y2,B.all_pos_y2);
% pos_t_all1=AppendTwoVectors(A.all_pos_t1,B.all_pos_t1);
% pos_t_all2=AppendTwoVectors(A.all_pos_t2,B.all_pos_t2);
% 
% PF_radii1=AppendTwoVectors(A.PF_radii1,B.PF_radii1);
% PF_radii2=AppendTwoVectors(A.PF_radii2,B.PF_radii2);
% gaussian_mat1=AppendTwoVectors(A.gaussian_mat1,B.gaussian_mat1);
% gaussian_mat2=AppendTwoVectors(A.gaussian_mat2,B.gaussian_mat2);
% 
% save('corr coef results of same arenas COMBINED', ...
%     'all_rates_1', 'all_rates_2',...
%     'all_zm_1', 'all_zm_2',...
%     'all_rm_1', 'all_rm_2',...
%     'all_corrs', 'all_corrs2', ...
%     'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',...
%    'pos_x_all1','pos_x_all2','pos_y_all1','pos_y_all2','pos_t_all1','pos_t_all2',...
%    'PF_radii1','PF_radii2','gaussian_mat1','gaussian_mat2');

load('simulated results rescaling arenas 300x per cell xenia.mat')
means1=nanmean(all_stability_corrs);
clearvars -except means1

load('simulated results rescaling arenas 200x per cell ohad.mat')
means2=nanmean(all_stability_corrs);
clearvars -except means1 means2 

load('simulated results rescaling arenas 300x per cell gilad.mat')
means3=nanmean(all_stability_corrs);
clearvars -except means1 means2 means3

load('simulated results rescaling arenas 200x per cell sophie.mat')
means4=nanmean(all_stability_corrs);
clearvars -except means1 means2 means3 means4 

%load('simulated fano factor means.mat')

combined_means=AppendTwoVectors(means1,means2) ;
combined_means=AppendTwoVectors(combined_means,means3) ;
combined_means=AppendTwoVectors(combined_means,means4) ;

save('simulated rescaling arena stability corr means','combined_means')

disp('')

function c=AppendTwoVectors(a,b)

%b=b(18:end);

if iscell(a)
    c=[a b];
else
    c(1:length(a))=a;
    c(end+1:end+length(b))=b;
end
