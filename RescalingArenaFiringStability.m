function RescalingArenaFiringStability

dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')

load('rescaling arenas data info.mat')
count=1;
grid_num=0;
MD_num=0;
stretch_num=0;
fields_num=0;
[all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2,all_rm_1, all_rm_2, all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,count,grid_num,MD_num,stretch_num,fields_num,...
    all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2]=...
    RescalingArenaStability(max_indices_all,peak_rates_all,zone_mats_all,...
    rate_mats_all,gridness_all, MD_scores_all,count,grid_num,MD_num,stretch_num,...
    fields_num,pos_t_all,pos_x_all,pos_y_all,autocorrs_all);

save('corr coef results of rescaled arenas', ...
    'all_rates_1', 'all_rates_2',...
    'all_zm_1', 'all_zm_2',...
    'all_rm_1', 'all_rm_2',...
    'all_corrs', 'all_corrs2', ...
    'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',...
    'all_pos_x1','all_pos_x2','all_pos_y1','all_pos_y2','all_pos_t1','all_pos_t2',...
    'gaussian_mat1','gaussian_mat2','PF_radii1','PF_radii2');

% for rescaling arenas, first and last always same

function [all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2, all_rm_1, all_rm_2,all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,count,grid_num,MD_num,stretch_num,fields_num,...
    all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2]=...
    RescalingArenaStability(max_indices_all,peak_rates_all,zone_mats_all,...
    rate_mats_all,gridness_all, MD_scores_all,count,grid_num,MD_num, stretch_num,...
    fields_num,pos_t_all,pos_x_all,pos_y_all,autocorrs_all)

for i=1:length(zone_mats_all)
    
    %properties of the cell
    max_indices= max_indices_all{i};
    peak_rates_all_arenas= peak_rates_all{i};
    zone_mats=zone_mats_all{i};
    rate_mats=rate_mats_all{i};
    gridness_scores=gridness_all{i};
    MDs=MD_scores_all{i};
    
    pos_x_aa=pos_x_all{i};
    pos_y_aa=pos_y_all{i};
    pos_t_aa=pos_t_all{i};
    
    %for PF radius
    autocorrs_aa=autocorrs_all{i};
    
    arena_len= length(max_indices)-1;
    
    for h=1:arena_len-1
        
        for k=h+1:arena_len
            
            ind1=h;
            ind2=k;
            
            grid1=gridness_scores(ind1);
            grid2=gridness_scores(ind2);
            MD1=MDs(ind1);
            MD2=MDs(ind2);
            
            if iscell(grid1)
                grid1=cell2mat(grid1);
                grid2=cell2mat(grid2);
                MD1=cell2mat(MD1);
                MD2=cell2mat(MD2);
            end
            
            if grid1> 0.3 || grid2>0.3
                
                grid_num=grid_num+1;
                
                if MD1<0.25 || MD2< 0.25
                    
                    MD_num=MD_num+1;
                    
                    max_inds_1=max_indices{ind1};
                    max_inds_2=max_indices{ind2};
                    
                     len_diff= length(max_inds_1)/length(max_inds_2);
                    if abs(1- (1/len_diff)) < 0.3 
                     
                        stretch_num=stretch_num+1;
                        
                    rates_1=peak_rates_all_arenas{ind1};
                    rates_2=peak_rates_all_arenas{ind2};
                    zm_1=zone_mats{ind1};
                    zm_2=zone_mats{ind2};
                    rm_1=rate_mats{ind1};
                    rm_2=rate_mats{ind2};
                    
                    size_1=size(zm_1);
                    size_2=size(zm_2);
                    new_size= [max(size_1(1),size_2(1)) max(size_1(2),size_2(2))];
                    max_inds_1(:,1)=max_inds_1(:,1)*(new_size(1)/size_1(1));
                    max_inds_1(:,2)=max_inds_1(:,2)*(new_size(2)/size_1(2));
                    max_inds_1=ceil(max_inds_1);
                    
                    [zm_1,zm_2]=ArenaSameSize(zm_1,zm_2);

                        if length(max_inds_1) >= 3
                            
                            fields_num=fields_num+1;
                            
                            % find peak rates of second session using 1st session inds
                            peak_rates_2= nan(1,length(max_inds_1));
                            for cen= 1:length(max_inds_1);
                                peak_rates_2(cen)= zm_2(max_inds_1(cen,1), max_inds_1(cen,2));
                            end
                            
                            peak_rates_1=rates_1;
                            
                            % find different correlations
                            corr_firing= corrcoef(peak_rates_1,peak_rates_2);
                            all_corrs(count)=corr_firing(2);
                            
                            prb_orig=peak_rates_1;
                            pre_orig=peak_rates_2;
                            
                            peak_rates_1(prb_orig==0 | pre_orig==0)= [];
                            peak_rates_2(prb_orig==0 | pre_orig==0)= [];
                            
                            if length(peak_rates_1) >2
                                corr_two= corrcoef(peak_rates_1,peak_rates_2);
                                all_corrs2(count)=corr_two(2);
                            else
                                all_corrs2(count)=nan;
                            end
                            
                            % twoD_zone_corrs(count)=corr2(zm_1,zm_2);

                            all_rates_1{count}=peak_rates_1;
                            all_rates_2{count}=peak_rates_2;
                            orig_rates_1{count}=rates_1;
                            orig_rates_2{count}=rates_2;
                            all_zm_1{count}=zm_1;
                            all_zm_2{count}=zm_2;
                            all_max_inds_1{count}=max_inds_1;
                            all_max_inds_2{count}=max_inds_2;
                            all_rm_1{count}=rm_1;
                            all_rm_2{count}=rm_2;
                            
                            %for simulations:
                    all_pos_x1{count}=pos_x_aa{ind1};
                    all_pos_x2{count}=pos_x_aa{ind2};
                    all_pos_y1{count}=pos_y_aa{ind1};
                    all_pos_y2{count}=pos_y_aa{ind2};
                    all_pos_t1{count}=pos_t_aa{ind1};
                    all_pos_t2{count}=pos_t_aa{ind2};

                    auto_max_inds1=FindAutoMaxInds(autocorrs_aa{ind1});
                    PF_rad1= findPlaceFieldRadius(autocorrs_aa{ind1}, auto_max_inds1)  ;
                    gaussian_mat1{count}=createGaussianMat(zm_1, PF_rad1,max_inds_1, mean(rates_1));
                    
                     auto_max_inds2=FindAutoMaxInds(autocorrs_aa{ind2});
                    PF_rad2= findPlaceFieldRadius(autocorrs_aa{ind2}, auto_max_inds2)  ;
                    gaussian_mat2{count}=createGaussianMat(zm_2, PF_rad2,max_inds_2, mean(rates_2));
                    
                    PF_radii1(count)=PF_rad1;
                    PF_radii2(count)=PF_rad2;
%                             figure;
%                             subplot(1,2,1)
%                             imagesc(zm_1)
%                             subplot(1,2,2)
%                             imagesc(zm_2)
                            
                            count=count+1;
                            
                            close all
                        end
                    end
                end
            end
        end
    end
end




