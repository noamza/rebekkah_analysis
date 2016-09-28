function hyperfield_remapping_analysis_UPDATED

dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info');

load('remapping data info.mat')

PF_thresh=3; % at least 3 fields
xcorr_thresh= 0;
phase_thresh= 0.2;

same_count=1;
remap_count=1;
no_remap_count=1;
all_count=1;

%find hyperfield distances between same context arena pairs
for h= 1:length(zone_mats_all)
    
    % take two hyperfield indices, find distance between them
    % find minimum distance, subtract dist from this amount to normalize
    
    % trial= inds(h);
    trial=h;
    
    % open all arena parameters
    hyper_inds_all_arenas= norm_max_index_all{trial};
    max_inds_all_arenas= max_indices_all{trial};
    rates_all_arenas= peak_rates_all{trial};
    zone_mats_all_arenas=zone_mats_all{trial};
    stability_inds= same_arena_inds{trial};
    rate_maps_all_arenas= rate_mats_all{trial};
    MD_scores_all_arenas= MD_scores_all{trial};
    gridness_all_arenas=gridness_scores_all{trial};
    autocorrs_all_arenas= autocorrs_all{trial};
    
    % check if same context trials is stable enough
    rate_map_1= rate_maps_all_arenas{stability_inds(1)};
    rate_map_2= rate_maps_all_arenas{stability_inds(2)};
    zone_mat_1= zone_mats_all_arenas{stability_inds(1)};
    zone_mat_2= zone_mats_all_arenas{stability_inds(2)};
    
    [zone_mat_1, zone_mat_2]= ...
        ArenaSameSize(zone_mat_1, zone_mat_2);
    
     [rate_map_1, rate_map_2]= ...
        ArenaSameSize(rate_map_1, rate_map_2);
    
    % convert zone mat to 0s and 1s
    zone_mat_1(zone_mat_1~=0)=1;
    zone_mat_2(zone_mat_2~=0)=1;
    
    xcorr= Cross_Correlation(zone_mat_1,zone_mat_2);
    R_outer=FindROuter(xcorr);
    gridness3= GridnessRadiusXcorr(xcorr,R_outer);
    
    if gridness3 > xcorr_thresh % if same context stable enough
        %
        %         acorr1= Cross_Correlation(rate_map_1,rate_map_1);
        %         spacing1= findPlaceFieldRadius(acorr1) * (10/7);
        %         acorr2= Cross_Correlation(rate_map_2,rate_map_2);
        %         spacing2= findPlaceFieldRadius(acorr2) * (10/7);
        %
        %         [phase_shift]= FindPhaseShift(xcorr, [spacing1 spacing2]);
        %
        %         if phase_shift < phase_thresh % if same context trial stable enough
        
        len=5;
        
        for a= 1:len-1
            for b= a+1:len
                
                rate_mat_1= rate_maps_all_arenas{a};
                rate_mat_2= rate_maps_all_arenas{b};
                
                max_inds_1= max_inds_all_arenas{a};
                max_inds_2= max_inds_all_arenas{b};
                
                rates_1=rates_all_arenas{a};
                rates_2=rates_all_arenas{b};
                
                zone_mat_1= zone_mats_all_arenas{a};
                zone_mat_2= zone_mats_all_arenas{b};
                
                hyper_ind_1= hyper_inds_all_arenas(a,:);
                hyper_ind_2= hyper_inds_all_arenas(b,:);
                
                [lenn_1,~] =size(max_inds_1);
                [lenn_2,~]= size(max_inds_2);
                
                %normalize max inds to arena size
                size_2= size(rate_mat_2);
                norm_max_inds_2=nan(lenn_2,2);
                norm_max_inds_2(:,1)= max_inds_2(:,1)/size_2(1);
                norm_max_inds_2(:,2)= max_inds_2(:,2)/size_2(2);
                
                MD_1= MD_scores_all_arenas{a};
                MD_2= MD_scores_all_arenas{b};
                
                grid1=gridness_all_arenas{a};
                grid2=gridness_all_arenas{b};
                
                acorr1= autocorrs_all_arenas{a};
                acorr2= autocorrs_all_arenas{b};
                
                % if enough number of fields in both trials and all other
                % parameters
                if lenn_1 >= PF_thresh && lenn_2 >= PF_thresh &&...
                        MD_1 <= 0.25 && MD_2 <= 0.25 &&...
                        grid1>=0.3 && grid2>=0.3
                    
                    dist=FindDistTwoHypers(hyper_ind_1,hyper_ind_2,norm_max_inds_2);
                    [gc_corr, orig, next]= FindFiringProfileGC(max_inds_1, max_inds_2, rates_1, rates_2);
                    
                    [zone_mat_1, zone_mat_2]= ...
                        ArenaSameSize(zone_mat_1, zone_mat_2);
                    
                    % need only for images:
                    %zone_mat_1_orig= zone_mat_1;
                    %zone_mat_2_orig= zone_mat_2;
                    
                    % convert zone mat to matrix of 0s and 1s
                    zone_mat_1(zone_mat_1~=0)=1;
                    zone_mat_2(zone_mat_2~=0)=1;
                    
                    % find phase shift
                    xcorr= Cross_Correlation(rate_mat_1,rate_mat_2);
                    xcorr_rate=xcorr;
                    
                    spacing1= findPlaceFieldRadius(acorr1) * (10/7);
                    spacing2= findPlaceFieldRadius(acorr2) * (10/7);
                    
                    [phase_shift]= FindPhaseShift(xcorr, [spacing1 spacing2]);
                    
                    xcorr= Cross_Correlation(zone_mat_1,zone_mat_2);
                    R_outer=FindROuter(xcorr);
                    [gridness3] = GridnessRadiusXcorr(xcorr,R_outer);
                    
                    all_dists(all_count)=dist;
                    all_corrcoefs(all_count)=gc_corr;
                    phase_shifts(all_count)= phase_shift;
                    xcorr_scores(all_count)=gridness3;
                    all_count=all_count+1;
                    
                    %divide by category
                    if isequal([a b], stability_inds)  % if same context trial pair
                        same_dists(same_count)=dist;
                        same_corrs(same_count)=gc_corr;
                        same_count=1;
                        
                    else  % if not same context pair
                        
                        if phase_shift < phase_thresh & gridness3 > xcorr_thresh %didnt remap
                            
                            no_remap_dists(no_remap_count)= dist;
                            no_remap_corrs(no_remap_count)=gc_corr;
                            no_remap_count=no_remap_count+1;
                            
                        elseif phase_shift >= phase_thresh | gridness3 <= xcorr_thresh
                            
                            remap_dists(remap_count)= dist;
                            remap_corrs(remap_count)=gc_corr;
                            remap_count=remap_count+1;
                            
                        end
                        
                        
                        disp('');
                        
                    end %of if same context
                end %of if both arenas pass all criteria                
            end % of arena 2
        end %of arena 1
        disp('');
    end %end of if same context phase stable
    %  end % end of if same context xcorr gridness stable 
end % of all cells

a=same_corrs;
b=no_remap_corrs;
c=remap_corrs;

fig= PlotMeanErrorbars(a,b,c);

a=same_dists;
b=no_remap_dists;
c=remap_dists;

fig2= PlotMeanErrorbars(a,b,c);


    



% bar graph of hyperfield location/category
nonremap_per= sum(b>0.4)/ length(b);
remap_per= sum(c>0.4)/ length(c);

nonremap_o= sum(b<=0.4)/ length(b);
remap_o= sum(c<=0.4)/ length(c);

y= [nonremap_per nonremap_o; remap_per remap_o];

figure; bar(y, 'stacked');
strng= sprintf('n=%d', sum(~isnan(a)));
text(1,1.2,strng);
strng= sprintf('n=%d', length(b));
text(2,1.2,strng);

%scatterplot of phase shift vs firing corr
figure; scatter(norm_shifts,all_corrcoefs)

figure; hist(norm_shifts);

dist_means= [mean(b) mean(c)];
stderr_a= FindStdErr(a);
stderr_b= FindStdErr(b);
stderr_c= FindStdErr(c);
dist_stderrs= [stderr_a stderr_b  stderr_c];

a= same_corrs;
b= no_remap_corrs;
c= remap_corrs;

corr_means= [nanmean(a) nanmean(b) nanmean(c)];
stderr_a= FindStdErr(a);
stderr_b= FindStdErr(b);
stderr_c= FindStdErr(c);
corr_stderrs= [stderr_a stderr_b  stderr_c];

figure;
errorbar(dist_means, dist_stderrs, 'x');
ylabel('Distance')
figure;
errorbar(corr_means, corr_stderrs, 'x');
ylabel('Corr Coef')

save('dists and corrs 2 2 3 phase', 'same_dists', 'no_remap_dists', 'remap_dists',...
    'same_corrs', 'no_remap_corrs', 'remap_corrs',...
    'dist_means', 'dist_stderrs', 'corr_means', 'corr_stderrs', ...
    'all_dists', 'all_corrcoefs', 'phase_shift', 'norm_shifts', 'all_2D_corrs')

disp('');

% bin data to see corr vs dist scatter

a= norm_shifts;

lims= [-0.4:0.2:1];

bins=nan(1,length(a));
for h=1:length(a)
    for k=1:length(lims)-1
        inds=(a >= lims(k) & a < lims(k+1));
        bins(inds)= k;
    end
end

b=all_corrcoefs;

for h=1:length(lims)-1
    x= b(bins==h);
    means(h)= nanmean(x);
    stderrs(h)= FindStdErr(x);
end

figure; errorbar(means, stderrs, 'ro-');

sum(isnan(all_corrcoefs))

save('binned data', 'means', 'stderrs', 'lims')

disp('')

remap= all_dists(xcorr_scores <= 0.1 | phase_shifts >= 0.1);
noremap= all_dists(xcorr_scores > 0.1 & phase_shifts <0.1);

remap_2corr= all_dists(all_2D_corrs <= 0.2);
noremap_2corr= all_dists(all_2D_corrs > 0.2);

a= remap_dists;
b= no_remap_dists;
xcorr_sums= [sum(a==0)/length(a) sum(a~=0)/length(a);
    sum(b==0)/length(b) sum(b~=0)/length(b)];
figure;bar(xcorr_sums, 'stacked')

a= remap_phase;
b= noremap_phase;
sums= [sum(a==0)/length(a) sum(a~=0)/length(a);
    sum(b==0)/length(b) sum(b~=0)/length(b)];
figure;bar(sums, 'stacked')

stderr1= FindStdErr(a);
stderr2= FindStdErr(b);
figure; errorbar([mean(a) mean(b)], [stderr1 stderr2], 'o')



function stderr= FindStdErr(x)

stderr= nanstd(x)/ sqrt((sum(~isnan(x))));

function fig= PlotMeanErrorbars(a,b,c)

    all_means= [mean(a) mean(b) mean(c)];
    stderr1=FindStdErr(a);
    stderr2=FindStdErr(b);
    stderr3=FindStdErr(c);
    all_stderrs= [stderr1 stderr2 stderr3];
    
    fig=figure; errorbar(all_means, all_stderrs)

