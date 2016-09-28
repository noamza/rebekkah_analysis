function SameArenaFiringStability

dbstop if error

%SameArenaStability();

% cd('C:\Noam\Output\rebekkah\');
% load('remapping data info noam mids');
% %load('remapping data info.mat')
% cd('C:\Noam\Dropbox\GitTechnion\rebekkah');

%function SameArenaStability
load('C:\Noam\Output\rebekkah\remapping data info noam mids');
count=1;
grid_num=0;
MD_num=0;
nonremap_num=0;
mi = 0; md = 0; gs = 0; va = 0; rm = 0;
peak_corrs=[]; all_corrs = [];
prnt = false;

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
    i
    %for PF radius
    autocorrs_aa=autocorrs_all{i}; h = [];
    h = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
    %figure;
    set(gca,'LooseInset', get(gca,'TightInset'));
    colormap jet;
    titl = sprintf('Cell %d', filenames{i}); %num2str(gridness_scores));
    set(suptitle(titl),'Interpreter', 'none');
    %valid_all{i}
    prnt = false;
    for j = 1:length(pos_x_aa) - 1;
        for k = j+1:3; %for k = j+1:length(pos_x_aa);
        %k = j + 1;
        ind1 = j; ind2 = k;
        if valid_all{i}{ind1} && valid_all{i}{ind2} 
            arena_type = arena_types_all{i};
            arena1 = arena_type{ind1};
            arena2 = arena_type{ind2};            
            grid1=gridness_scores(ind1);
            grid2=gridness_scores(ind2);
            MD1=MDs(ind1);
            MD2=MDs(ind2);
            %if grid1> 0.3 || grid2>0.3 || true %NOAM
            if grid1> 0.5 && grid1 - grid2 > 0.2 && ind1 == 1;
                grid_num=grid_num+1;
                if MD1<0.15 || MD2< 0.15
                    MD_num=MD_num+1;
                    max_coords_1=max_indices{ind1};
                    max_coords_2=max_indices{ind2};
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
                    zm101=zm_1;zm101(zm101~=0)=1;zm102=zm_2;zm102(zm102~=0)=1;
                    remap_level=corr2(zm101,zm102);
                    remap_level = 999999; %NOAM
                    if remap_level > 0.3 %remap
                        nonremap_num=nonremap_num+1;
                        if length(max_coords_1) >= 3 && length(max_coords_2) >= 3
                            prnt = true;                  %%% PRINT? %%%
                            % find peak rates of second session using 1st session inds
                            %sort by largest peak
                            [peak_rates_1_ordered,order] = sort(rates_1, 'descend'); %order by peak rate of 1
                            %indices in order of peak rate of 1
                            max_coords_1_ordered = max_coords_1(order,:);
                            %peak_rates_2 = findPeakRates(max_coords_1, zm_2); %in order
                            [indices, pairs, distances] = NoamFindClosestPeaks(max_coords_1_ordered, max_coords_2);
                            %peak_rates_2 = findPeakRates(max_coords_1,%zm_2); 
                            %sort everything by firing rates of 1
                            %works because original indexing order is based on rates_1 ordered
                            order = []; [indices, order] = sortrows(indices,1); %sorted by rates_1 
                            pairs = pairs(order,:);
                            distances = distances(order);
                            peak_rates_1_ordered = peak_rates_1_ordered(indices(:,1)); %only uses rates used in pairs
                            peak_rates_2_ordered_by_r1 = rates_2(indices(:,2)); %only uses rates used in pairs, ordered by rates 1
                            %figure; scatter(1:5,peak_rates_1_ordered); hold on; scatter(1:5,peak_rates_2_ordered_by_r1);
                            %EACH DIRECTION
                            % find different correlations
                            corr_firing = corrcoef(peak_rates_1_ordered, peak_rates_2_ordered_by_r1);
                            if( k == 2)
                                peak_corrs(end + 1) = corr_firing(2);
                            end
                            all_corrs(end + 1)=corr_firing(2); % (count);
                            %titl = sprintf('cell %d: %s vs. %s', filenames{i}, arena1, arena2);
                            %%%%NOAM%%%%
                            cols = 6; %number of cols
                            rows = 2;
                            %rate mat 1
                            subplot(rows,cols,(k-2)*cols+1);
                            colorbar;
                            imagesc(rate_mats{ind1});
                            set(gca,'YDir','normal');
                            title(sprintf('%s (gs %.1f)',arena1, gridness_scores(ind1)));
                            axis equal; axis off;
                            %rate mat 2
                            subplot(rows,cols,(k-2)*cols+2);
                            imagesc(rate_mats{ind2});
                            set(gca,'YDir','normal');
                            title(sprintf('%s (gs %.1f)',arena2,gridness_scores(ind2)));
                            axis equal; axis off;
                            %zone 1
                            subplot(rows,cols,(k-2)*cols+3);
                            imagesc(zm_1);
                            set(gca,'YDir','normal');
                            hold on ; colormap jet
                            %text(max_inds_1(:,2),max_inds_1(:,1), num2str((1:length(max_inds_1))'), ...
                            %text(pairs(:,2),pairs(:,1), num2str(peak_rates_1'), ...
                            text(pairs(:,2),pairs(:,1), num2str((1:length(distances))'), ...
                                'HorizontalAlignment','center','VerticalAlignment','middle','Color','w','FontSize',12,'FontWeight','bold');
                            title(sprintf('%s (%.1fHz)',arena1, max(peak_rates_1_ordered)));
                            %subtitle(sprintf('%s %.1fHz','max',max(peak_rates_1)));
                            axis equal;  axis off;
                            %zone mat 2
                            subplot(rows,cols,(k-2)*cols+4);
                            imagesc(zm_2);
                            set(gca,'YDir','normal');
                            hold on; colormap jet
                            title(sprintf('%s vs %s (%.1fHz)',arena1, arena2, max(peak_rates_2_ordered_by_r1)));
                            %text(max_inds_1(:,2),max_inds_1(:,1), num2str((1:length(max_inds_1))'), ...
                            text(pairs(:,4),pairs(:,3), num2str((1:length(distances))'), ...
                                'HorizontalAlignment','center','VerticalAlignment','middle','Color','w','FontSize',12,'FontWeight','bold');
                            axis equal;  axis off;
                            
                            %zone compare
                            subplot(rows,cols,(k-2)*cols + 5); hold on
                            scatter(pairs(:,2),pairs(:,1),120,linspace(0,1,length(pairs(:,1))),'filled');
                            scatter(pairs(:,4),pairs(:,3),120,linspace(0,1,length(pairs(:,1))),'filled');
                            text(pairs(:,2),pairs(:,1), '1', 'HorizontalAlignment','right','VerticalAlignment','bottom','Color','black','FontSize',12,'FontWeight','bold');
                            text(pairs(:,4),pairs(:,3), '2', 'HorizontalAlignment','right','VerticalAlignment','bottom','Color','black','FontSize',12,'FontWeight','bold');
                            title(sprintf('%s','Zone Pairs'));
                            axis equal tight;
                            
                            %scatter firint rates
                            subplot(rows,cols,(k-2)*cols+6);
                            scatter(peak_rates_1_ordered, peak_rates_2_ordered_by_r1);
                            title(sprintf('peak corr: %.2f',corr_firing(2)));
                            axis equal square;
                            xlabel('firing rate 1');
                            ylabel('firing rate 2');
                           
 
                            ax = findobj(gcf,'Type','Axes');
                            for z=1:length(ax)
                                %set(ax(z),'FontSize',9);
                                %axis(ax(z),'equal')
                                %axis(ax(z),'off')
                                %axis off; axis equal;
                                %title(ax(z),{'Very Nice'})
                            end
                            disp('');
                            %%%%%%%%%
                            
                            % FROM THE MIDDLE
                        else
                            mi = mi + 1;
                            %disp('not analyzed because max inds < 3');
                        end
                    else
                        rm = rm + 1;
                        %disp('not analyzed because remapped');
                    end
                else
                     md = md + 1; 
                    %disp('not analyzed because md small');
                end
            else
                gs = gs + 1;
                %disp('not analyzed because invalid');
            end
        else
             va = va + 1;
             %disp('not analyzed because invalid');
        end
      end
    end
    if prnt
        debug = '';
        sdir = 'C:\Noam\Output\rebekkah\plots\';
        filename = sprintf('%s%s%s.png',sdir, debug, titl);
        %titl = sprintf('cells%d_pg%dof%d\n', length(g), page, ceil(length(g)/r));
        disp(filename);
        h.PaperPositionMode = 'auto'; %fig %fig = gcf;
        print(filename, '-dpng','-r0');
    end
    close(h);
end
fprintf('max peaks < 3(%d)\nmoving direction < .15(%d)\nremap level < 0.3(%d)\ngridscore < 0.3(%d)\ninvalid from preprocessing(%d)\n',... 
    mi,md,rm,gs,va); 

h = figure();
hist(peak_corrs,10);
title('Peak correlations before-muscimol grid1-grid2 > 0.2');
sdir = 'C:\Noam\Output\rebekkah\plots\';
filename = sprintf('%s%s',sdir,'hist corrs before mid1.png');
disp(filename);
h.PaperPositionMode = 'auto'; %fig %fig = gcf;
print(filename, '-dpng','-r0');

h = figure();
hist(all_corrs,10);
title('All peak correlations before-all grid1-grid2 > 0.2');
sdir = 'C:\Noam\Output\rebekkah\plots\';
filename = sprintf('%s%s',sdir,'all hist corrs.png');
disp(filename);
h.PaperPositionMode = 'auto'; %fig %fig = gcf;
print(filename, '-dpng','-r0');

end

%end
%CUT

%{ from the middle

                            %
                            %                             prb_orig = peak_rates_1;
                            %                             pre_orig = peak_rates_2;
                            %
                            %                             peak_rates_1(prb_orig==0 | pre_orig==0)= [];
                            %                             peak_rates_2(prb_orig==0 | pre_orig==0)= [];
                            %
                            %                             if length(peak_rates_1) >2
                            %                                 corr_two= corrcoef(peak_rates_1,peak_rates_2);
                            %                                 all_corrs2(count)=corr_two(2);
                            %                             else
                            %                                 all_corrs2(count)=nan;
                            %                             end
                            %
                            %                             % twoD_zone_corrs(count)=corr2(zm_1,zm_2);
                            %
                            %                             all_rates_1{count}=peak_rates_1;
                            %                             all_rates_2{count}=peak_rates_2;
                            %                             orig_rates_1{count}=rates_1;
                            %                             orig_rates_2{count}=rates_2;
                            %                             all_zm_1{count}=zm_1;
                            %                             all_zm_2{count}=zm_2;
                            %                             all_max_inds_1{count}=max_inds_1;
                            %                             all_max_inds_2{count}=max_inds_2;
                            %                             all_rm_1{count}=rm_1;
                            %                             all_rm_2{count}=rm_2;
                            %
                            %                             %for simulations:
                            %                             all_pos_x1{count}=pos_x_aa{ind1};
                            %                             all_pos_x2{count}=pos_x_aa{ind2};
                            %                             all_pos_y1{count}=pos_y_aa{ind1};
                            %                             all_pos_y2{count}=pos_y_aa{ind2};
                            %                             all_pos_t1{count}=pos_t_aa{ind1};
                            %                             all_pos_t2{count}=pos_t_aa{ind2};
                            %
                            %                             auto_max_inds1=FindAutoMaxInds(autocorrs_aa{ind1});
                            %                             PF_rad1= findPlaceFieldRadius(autocorrs_aa{ind1}, auto_max_inds1)  ;
                            %                             gaussian_mat1{count}=CreateGaussianMat(zm_1, PF_rad1,max_inds_1, mean(rates_1));
                            %
                            %                             auto_max_inds2=FindAutoMaxInds(autocorrs_aa{ind2});
                            %                             PF_rad2= findPlaceFieldRadius(autocorrs_aa{ind2}, auto_max_inds2)  ;
                            %                             gaussian_mat2{count}=CreateGaussianMat(zm_2, PF_rad2,max_inds_2, mean(rates_2));
                            %
                            %                             PF_radii1(count)=PF_rad1;
                            %                             PF_radii2(count)=PF_rad2;
                            %                             %                     figure;
                            %                             %                     subplot(1,2,1)
                            %                             %                     imagesc(zm_1)
                            %                             %                     subplot(1,2,2)
                            %                             %                     imagesc(zm_2)
                            %                             count=count+1;
                            %
                            %                             close all

%}


%{ 
from the top


count=1;
grid_num=0;
MD_num=0;
nonremap_num=0;

% cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
% load('rescaling arenas data info.mat')
% [all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
%     all_zm_1, all_zm_2,all_rm_1, all_rm_2, all_max_inds_1, all_max_inds_2,...
%     all_corrs, all_corrs2,...
%     all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
%     gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2,...
%     count,grid_num,MD_num,nonremap_num]=...
%     SameArenaStability(pos_x_all,pos_y_all,pos_t_all,autocorrs_all,...
%     max_indices_all,peak_rates_all,zone_mats_all,...
%     rate_mats_all,gridness_all, MD_scores_all,count,grid_num,MD_num,nonremap_num);
%
% save('same arenas corr coef results of rescaling data', ...
%     'all_rates_1', 'all_rates_2',...
%     'all_zm_1', 'all_zm_2',...
%      'all_rm_1', 'all_rm_2',...
%     'all_corrs', 'all_corrs2', ...
%     'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',   'all_pos_x1', 'all_pos_x2','all_pos_y1','all_pos_y2', 'all_pos_t1','all_pos_t2',...
%     'gaussian_mat1','gaussian_mat2','PF_radii1','PF_radii2');
%
% clearvars -except count grid_num MD_num nonremap_num
% load('remapping data info.mat')

cd('C:\Noam\Output\rebekkah\');
load('remapping data info noam mids');
%load('remapping data info.mat')
cd('C:\Noam\Dropbox\GitTechnion\rebekkah');    
    
    
same_arena_inds=cell(1,length(arena_types_all));
% for c = 1:length(arena_types_all)
%
%     arena_type= arena_types_all{c};
%
%     len = length(arena_type); %5
%     for a= 1:len-1
%         for b=a+1:len
%             A= arena_type{a};
%             B= arena_type{b};
%             %if A(1:2) == B(1:2)
%             if true %noam
%                 sa_inds=[a b];
%             end
%         end
%     end
%
%     same_arena_inds{c}= sa_inds;
% end


[all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2,all_rm_1, all_rm_2, all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,...
    all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2,...
    count,grid_num,MD_num,nonremap_num]=...
    SameArenaStability(pos_x_all,pos_y_all,pos_t_all,autocorrs_all,...
    max_indices_all, peak_rates_all,zone_mats_all,...
    rate_mats_all,gridness_all, MD_scores_all,count,grid_num,MD_num,nonremap_num,...
    same_arena_inds, filenames, arena_types_all, valid_all);

save('same arenas corr coef results', ...
    'all_rates_1', 'all_rates_2',...
    'all_zm_1', 'all_zm_2',...
    'all_rm_1', 'all_rm_2',...
    'all_corrs', 'all_corrs2', ...
    'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',...
    'all_max_inds_1','all_max_inds_2','orig_rates_1','orig_rates_2',   'all_pos_x1', 'all_pos_x2','all_pos_y1','all_pos_y2', 'all_pos_t1','all_pos_t2',...
    'gaussian_mat1','gaussian_mat2','PF_radii1','PF_radii2');

% for rescaling arenas, first and last always same

SameArenaStability

function [all_rates_1, all_rates_2, orig_rates_1, orig_rates_2, ...
    all_zm_1, all_zm_2, all_rm_1, all_rm_2,all_max_inds_1, all_max_inds_2,...
    all_corrs, all_corrs2,...
    all_pos_x1,all_pos_x2,all_pos_y1,all_pos_y2,all_pos_t1,all_pos_t2,...
    gaussian_mat1,gaussian_mat2,PF_radii1,PF_radii2,...
    count,grid_num,MD_num,nonremap_num]=...
    SameArenaStability(pos_x_all,pos_y_all,pos_t_all,autocorrs_all,...
    max_indices_all,peak_rates_all,zone_mats_all,...
    rate_mats_all, gridness_all, MD_scores_all,count,grid_num,MD_num, nonremap_num,...
    same_arena_inds, filenames, arena_types_all, valid_all)
    
    
%}

