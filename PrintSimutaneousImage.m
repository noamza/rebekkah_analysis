function PrintSimutaneousImage(max_indices_sm,letter, inds1, rate_mats_all, step, com, rats_use,dates_use, sesh_use, peak_rates_sm_all)

for c=1:length(inds1)
    length(inds1)
    
    subplot(6,6,c+(step(c)))
    e=(inds1(c));
    imagesc(rate_mats_all{e})
    axis off; axis image; axis ij;
    if c==1
        title(sprintf('rat=%s, date=%s, session=%s', rats_use{e},dates_use{e},sesh_use{e}));
        text(-0.25,1.05,letter{1},'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
    end
    
    peak_rates=peak_rates_sm_all{e};
    max_ind=find(peak_rates==max(peak_rates));
    norm_ind= max_indices_sm{e};
    norm_ind=norm_ind(max_ind,:);
    hold on;
    plot(norm_ind(2),norm_ind(1),'ko', 'MarkerSize', 15, 'LineWidth', 2)
    
    
    rm_size= size(rate_mats_all{e});
    subplot(6,6,com)
    plot(norm_ind(2),norm_ind(1),'kx', 'LineWidth', 2)
    hold on;
    axis ij;
    axis square;
    xlim([0 rm_size(2)])
    ylim([0 rm_size(1)])
    %axis off;
    
end