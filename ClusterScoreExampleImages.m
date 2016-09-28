cd('\\192.114.21.198\Dori_Docs\users\rebekkah\results and info of analysis')
load('Cluster scores.mat')

cluster_1_inds= find(cluster_scores_35==1);
cluster_pt_6_inds= find(cluster_scores_35==0.6);
cluster_pt_4_inds= find(cluster_scores_35==0.4);
cluster_pt_3_inds= find(cluster_scores_35==0.3);
cluster_pt_1_inds= find(cluster_scores_35==0.1);
cluster_0_inds= find(cluster_scores_35==0);

cluster_1_inds(1)=[];
cluster_pt_3_inds(1:2)=[];

% cluster_inds(1:3)= cluster_1_inds(1:3);
% cluster_inds(4:5)= cluster_pt_6_inds(1:2);
% cluster_inds(6)= cluster_pt_4_inds(1);
% cluster_inds(7)= cluster_pt_3_inds(1);
% cluster_inds(8)= cluster_pt_1_inds(1);
% cluster_inds(9)= cluster_0_inds(1);

cluster_inds(1:2)= cluster_1_inds(1:2);
cluster_inds(3:4)= cluster_pt_6_inds(1:2);
%cluster_inds(6)= cluster_pt_4_inds(1);
cluster_inds(5)= cluster_pt_3_inds(1);
cluster_inds(6)= cluster_pt_1_inds(1);
cluster_inds(7)= cluster_0_inds(1);

load('rescaling arenas info.mat')

add=0;

letter= [{'A'} {'B'} {'C'} {'D'} {'E'} {'F'} {'G'} {'H'} {'I'} {'J'}];
scores= [{'1.0'} {'1.0'} {'0.6'} {'0.6'} {'0.3'} {'0.1'} {'0.0'}];

%plot example of methods
% fig=figure;
% n=7;
% m=2;
% 
% count=1;
% for h=1;
%     rate_mats= all_rate_mats{cluster_inds(h)};
%     norm_indices= all_norm_inds{cluster_inds(h)};
%     peak_rates_all_arenas= all_peak_rates{cluster_inds(h)};
%     max_indices_all_arenas= all_max_inds{cluster_inds(h)};
%     
%     max_index= nan(5,2);
%     for arena_count=1:5;
%         
%         peak_rates=peak_rates_all_arenas{arena_count};
%         norm_inds= norm_indices{arena_count};
%         max_ind= find(peak_rates==max(peak_rates));
%         max_index(arena_count,:)=norm_inds(max_ind, :);
%         max_indices= max_indices_all_arenas{arena_count};
%         hyperfield= max_indices(max_ind,:);
%         rate_mat=rate_mats{arena_count};
%         
%         size_rm= size(rate_mat);
%         
%         subplot(n,m,count)
%     %    plot(hyperfield(2),hyperfield(1), 'x', 'LineWidth', 12, 'MarkerSize', 3, 'color', 'k');
%         axis image; 
%         set(gca,'xtick',[], 'ytick',[]);
%         xlim([1 size_rm(2)]);
%         ylim([1 size_rm(1)]);
%         
%         subplot(n,m,count+1)
%         plot(max_index(arena_count,2),max_index(arena_count,1), 'x','LineWidth', 12, 'MarkerSize', 3, 'color', 'k');
%         axis image; 
%         set(gca,'xtick',[], 'ytick',[]);
%         xlim([0 1]);
%         ylim([0 1]);
%         
%         count=count+2;
%     end 
%     
%     subplot(n,m,count+3)
%     plot(max_index(:,2), max_index(:,1), 'x',  'LineWidth', 12, 'MarkerSize', 3, 'color', 'k');
%     set(gca,'xtick',[], 'ytick',[]);
%     axis image;
%     xlim([0 1]);
%     ylim([0 1]);
%   %  text(0.35,1,scores(h),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 15, 'Color', 'r')
%     
%     
% end
% 
% dist_bw_arrows= 0;
% for h=1:5
%     annotation(fig,'arrow',[0.47 0.56],...
%         [0.88+dist_bw_arrows 0.88+dist_bw_arrows]);
%     
%     dist_bw_arrows= dist_bw_arrows- 0.12;
% end
% 
%  annotation(fig,'arrow',[0.73 0.73],[0.3 0.24]);
% 
% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 24], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 24])

figure;
n= 7;
m= 6;
for h= 1:7
    rate_mats= all_rate_mats{cluster_inds(h)};
    norm_indices= all_norm_inds{cluster_inds(h)};
    peak_rates_all_arenas= all_peak_rates{cluster_inds(h)};
    max_indices_all_arenas= all_max_inds{cluster_inds(h)};
    
    max_index=nan(5,2);
    for arena_count=1:5
                
        peak_rates=peak_rates_all_arenas{arena_count};
        norm_inds= norm_indices{arena_count};
        max_ind= find(peak_rates==max(peak_rates));
        max_index(arena_count,:)=norm_inds(max_ind, :);
        max_indices= max_indices_all_arenas{arena_count};
        hyperfield= max_indices(max_ind,:);
        
        %plots rate mat wth hyperfield circled
        subplot(n,m,arena_count+add)
        %imagesc(rate_mats{arena_count}, 'AlphaData', 0.7); hold on;
        imagesc(rate_mats{arena_count}); hold on;
        axis equal; axis off;
        plot(hyperfield(2), hyperfield(1), 'o',  'LineWidth', 2, 'MarkerSize', 15, 'color', 'k');
        
        if arena_count==1
            load('file names rescaling data.mat')
            title(sprintf('%s', filenames{cluster_inds(h)}));
        end
        
        if arena_count==1
            text(-0.35,1.05,letter(h),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14, 'Rotation', 0)
        end
        
        
        %         subplot(n,m,arena_count+add+7)
        %         imagesc(rate_mats{arena_count}, 'AlphaData', 0.5);
        %         hold on; axis equal; axis off;
        %         plot(hyperfield(2), hyperfield(1), 'x',  'LineWidth', 8, 'MarkerSize', 3, 'color', 'k');
        
    end
    plot_inds=6 + add;
    
    %plots cluster score example 
    subplot(n,m,plot_inds)
    
 
    
    plot(max_index(:,2), max_index(:,1), 'x',  'LineWidth', 10, 'MarkerSize', 3, 'color', 'k');
    set(gca,'xtick',[], 'ytick',[]);
    axis square; axis ij;
    xlim([0 1]);
    ylim([0 1]);
    
 
    text(0.35,1,scores(h),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 15, 'Color', 'r', 'Rotation', 0)
    whitebg([0.8 0.8 0.8]);
    
    add=add+6;
end


set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 27.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 27.4])