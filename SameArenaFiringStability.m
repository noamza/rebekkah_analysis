function SameArenaFiringStability

dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')

load('rescaling arenas data info.mat')
count=1;
grid_num=0;
MD_num=0;
nonremap_num=0;
[all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2,all_rm_1, all_rm_2, all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,...
    all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2,...
    count,grid_num,MD_num,nonremap_num]=...
    SameArenaStability(pos_x_all,pos_y_all,pos_t_all,autocorrs_all,...
    max_indices_all,peak_rates_all,zone_mats_all,...
    rate_mats_all,gridness_all, MD_scores_all,count,grid_num,MD_num,nonremap_num);

save('same arenas corr coef results of rescaling data', ...
    'all_rates_1', 'all_rates_2',...
    'all_zm_1', 'all_zm_2',...
     'all_rm_1', 'all_rm_2',...
    'all_corrs', 'all_corrs2', ...
    'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',   'all_pos_x1', 'all_pos_x2','all_pos_y1','all_pos_y2', 'all_pos_t1','all_pos_t2',...
    'gaussian_mat1','gaussian_mat2','PF_radii1','PF_radii2');

clearvars -except count grid_num MD_num nonremap_num

load('remapping data info.mat')

same_arena_inds=cell(1,length(arena_types_all));
for c = 1:length(arena_types_all)
    
    arena_type= arena_types_all{c};
    
    len=5;
    for a= 1:len-1
        for b=a+1:len
            A= arena_type{a};
            B= arena_type{b};
            if A(1:2) == B(1:2)
                sa_inds=[a b];
            end
        end
    end
    
    same_arena_inds{c}= sa_inds;
end

[all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2,all_rm_1, all_rm_2, all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,...
   all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2,...
    count,grid_num,MD_num,nonremap_num]=...
    SameArenaStability(pos_x_all,pos_y_all,pos_t_all,autocorrs_all,...
    max_indices_all, peak_rates_all,zone_mats_all,...
    rate_mats_all,gridness_all, MD_scores_all,count,grid_num,MD_num,nonremap_num,...
    same_arena_inds);

save('same arenas corr coef results', ...
    'all_rates_1', 'all_rates_2',...
    'all_zm_1', 'all_zm_2',...
     'all_rm_1', 'all_rm_2',...
    'all_corrs', 'all_corrs2', ...
    'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',...
      'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',   'all_pos_x1', 'all_pos_x2','all_pos_y1','all_pos_y2', 'all_pos_t1','all_pos_t2',...
    'gaussian_mat1','gaussian_mat2','PF_radii1','PF_radii2');

% for rescaling arenas, first and last always same

function [all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2, all_rm_1, all_rm_2,all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,...
    all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2,...
    count,grid_num,MD_num,nonremap_num]=...
    SameArenaStability(pos_x_all,pos_y_all,pos_t_all,autocorrs_all,...
    max_indices_all,peak_rates_all,zone_mats_all,...
    rate_mats_all, gridness_all, MD_scores_all,count,grid_num,MD_num, nonremap_num,...
    same_arena_inds)

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
    
    if ~exist('same_arena_inds','var')
        len=length(max_indices);
        ind1=1;
        ind2=len;
    else
        sa_inds=same_arena_inds{i};
        ind1=sa_inds(1);
        ind2=sa_inds(2);
    end
    
    grid1=gridness_scores(ind1);
    grid2=gridness_scores(ind2);
    MD1=MDs(ind1);
    MD2=MDs(ind2);
    
    if grid1> 0.3 || grid2>0.3
        
        grid_num=grid_num+1;
        
        if MD1<0.15 || MD2< 0.15
            
            MD_num=MD_num+1;
            
            max_inds_1=max_indices{ind1};
            max_inds_2=max_indices{ind2};
            rates_1=peak_rates_all_arenas{ind1};
            rates_2=peak_rates_all_arenas{ind2};
            zm_1=zone_mats{ind1};
            zm_2=zone_mats{ind2};
            rm_1=rate_mats{ind1};
            rm_2=rate_mats{ind2};
            
            %             size_1=size(zm_1);
            %             size_2=size(zm_2);
            %             if ~isequal(size_1, size_2)
            %                 [zm_2] = StretchImage(zm_2, size_2, size_1);
            %             end
            
            [zm_1,zm_2]=ArenaSameSize(zm_1,zm_2);
            
            zm101=zm_1;
            zm101(zm101~=0)=1;
            zm102=zm_2;
            zm102(zm102~=0)=1;
            
            remap_level=corr2(zm101,zm102);
            
            if remap_level > 0.3
                
                nonremap_num=nonremap_num+1;
                
                if length(max_inds_1) >= 3
                    
                    % find peak rates of second session using 1st session inds
                    peak_rates_2=findPeakRates(max_inds_1,zm_2);
                    
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
                    gaussian_mat1{count}=CreateGaussianMat(zm_1, PF_rad1,max_inds_1, mean(rates_1));
                    
                     auto_max_inds2=FindAutoMaxInds(autocorrs_aa{ind2});
                    PF_rad2= findPlaceFieldRadius(autocorrs_aa{ind2}, auto_max_inds2)  ;
                    gaussian_mat2{count}=CreateGaussianMat(zm_2, PF_rad2,max_inds_2, mean(rates_2));
                    
                    PF_radii1(count)=PF_rad1;
                    PF_radii2(count)=PF_rad2;
                    %                     figure;
%                     subplot(1,2,1)
%                     imagesc(zm_1)
%                     subplot(1,2,2)
%                     imagesc(zm_2)
                    
                    count=count+1;
                    
                    close all
                end
            end
        end
    end
end




