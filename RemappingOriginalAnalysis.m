function RemappingOriginalAnalysis
%Date: 01 of June, 2016. Updated by Rebekkah.

dbstop if error

%for remapping data:
% parms.dir_load_data=('\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells with no pos data fixing');
%parms.dir_save_images=('\\192.114.21.198\Dori_Data\data\rebekkah\All rats\images original analysis');


% cd('\\192.114.21.198\Dori_Data\data\rebekkah\All rats\results and info of analysis');
% load('variability and border distances MD25 G3.mat', ...
%     'all_filenames', 'rate_mats_all', 'autocorrs_all', 'max_indices_all',...
%     'MD_scores_all', 'gridness_scores_all')

%for rescaling data:
% parms.dir_load_data=('\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells with no pos data fixing');
%

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info');
load('remapping data info.mat')
% load('Barry data info.mat')

if ~exist('filenames','var')
    filenames=all_filenames;
end

if ~exist('gridness_all','var')
    gridness_all=gridness_scores_all;
end

grid_num=0;
MD_num=0;
field_num=0;

c=1;

for i =1:length(rate_mats_all)
    
    posx_aa=pos_x_all{i};
    posy_aa=pos_y_all{i};
    post_aa=pos_t_all{i};
    spkx_aa=spk_x_all{i};
    spky_aa=spk_y_all{i};
    spkt_aa=spk_t_all{i};
    
    rate_mat_all_arenas= rate_mats_all{i};
    autocorrs_all_arenas= autocorrs_all{i};
    
    % MD_scores_all_arenas= MD_scores_all{i};
    max_indices_all_arenas=max_indices_all{i};
    gridness_all_arenas=gridness_all{i};
    MD_scores_all_arenas=MD_scores_all{i};
    filenames_all_arenas= filenames{i};
    
    PF_num=nan(1,length(max_indices_all_arenas));
    for h=1:length(max_indices_all_arenas)
        max_indices=max_indices_all_arenas{h};
        PF_num(h)= length(max_indices);
    end
    
    % check that all parameters pass criteria
    % med_MD_score= median(cell2mat(MD_scores_all_arenas)); % only for
    % remapping
    %[gridness_max, ind]=max(gridness_all_arenas); %for rescaling data
    
    %for remapping data:
    gridness_all_arenas=cell2mat(gridness_all_arenas);
    MD_scores_all_arenas=cell2mat(MD_scores_all_arenas);
    
    ind= find(gridness_all_arenas>0.3 & MD_scores_all_arenas<0.25 ...
        & PF_num >= 7);
    
  
    if ~isempty(ind)
        
         [len_ind,~]=size(ind);
    ind_lens(c)=len_ind;
        
        ind=ind(1);
        
        %         grid_num=grid_num+1;
        %
        %         if min_MD <=0.3
        %
        %             MD_num=MD_num+1;
        %
        %         if max_PF_num >= 7 %med_MD_score < 0.25 &&
        % cd(parms.dir_load_data)
        % filename= sprintf('%s%d.mat',filenames_all_arenas(1:end-5),ind);
        
        %         pos_x= (Cell.pos_data.x1);
        %         pos_y= (Cell.pos_data.y1);
        %
        %         pos_t= Cell.pos_data.t;
        %         spk_t= Cell.spike_data.ts;
        %
        %         spk_x=interp1(pos_t,pos_x,spk_t);
        %         spk_y=interp1(pos_t,pos_y,spk_t);
        
        field_num=field_num+1;
        
        pos_x= posx_aa{ind};
        pos_y= posy_aa{ind};
        pos_t= post_aa{ind};
        spk_x= spkx_aa{ind};
        spk_y= spky_aa{ind};
        spk_t= spkt_aa{ind};
        
        rate_mat=rate_mat_all_arenas{ind};
        autocorr= autocorrs_all_arenas{ind};
        
        [fanos(c), norm_max_index_ns(c,:), max_indices_i{c}, ...
            max_indices_sm_i{c}, norm_dist(c), rm_size(c,:), peak_rates_all_i{c}, ...
            peak_rates_all_sm{c}, norm_second_dist(c), ~,...
            ~, zone_mat_ns{c}]= ...
            findAllInfo(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t,rate_mat,autocorr);
        
        sorted_rates= sort(peak_rates_all_sm{c});
        max_over_means(c)= sorted_rates(end)/mean(sorted_rates(1:end-1));
        
        %         fig=figure; n=1; m=3;
        
        %         subplot(n,m,1)
        %         plot(pos_y,pos_x,'k'); hold on;
        %         plot(spk_y,spk_x,'r.'); axis off; axis equal;
        %         subplot(n,m,2)
        %         imagesc(rate_mat); axis off; axis square;
        %         title('smoothed rate mat')
        %         subplot(n,m,3);
        %         imagesc(zone_mat_ns{c}); axis off; axis square;
        %         title(sprintf('nonsmooth zone mat \n dist=%0.2f', norm_dist(c)));
        %
        %
        %         cd(parms.dir_save_images)
        %         saveas(fig,sprintf('%s.jpg', filename));
        %         cd(parms.dir_load_data)
        
        c=c+1;
        
        
        disp('')
    end
end



save('remapping data results', 'fanos', 'norm_dist', 'max_over_means', 'max_indices_i',...
    'rm_size', 'norm_second_dist', 'ind_lens')

