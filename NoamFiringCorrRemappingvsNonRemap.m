function FiringCorrRemappingvsNonRemap

binsize = 10;                            % B I N S I Z E
    
fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat',binsize);

dbstop if error

%cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
cd('C:\Noam\Output\rebekkah\');
load('remapping data info noam mids');
%load('remapping data info.mat')
cd('C:\Noam\Dropbox\GitTechnion\rebekkah');
same_arena_inds=cell(1,length(arena_types_all));
%finds indices of arenas with same 2 letter prefix WL or WV .. etc
for h=1:length(arena_types_all)
    types=arena_types_all{h};
    for a1=1:length(types)-1
        for a2=a1+1:length(types)
            type1=types{a1};
            type2=types{a2};
            %if type1(1:2)==type2(1:2)
            if true
                same_arena_inds{h}=[a1,a2];
            end
        end
    end
end


PF_thresh=3; % at least 3 fields
xcorr_thresh= 0;
phase_thresh= 0.2;

remap_count=1;
no_remap_count=1;
all_count=1;

no_remap_corrs = [] %NOAM
remap_corrs = []

for h= 1:length(zone_mats_all)

    % trial= inds(h);
    trial=h;

    % open all arena parameters
    %     hyper_inds_all_arenas= norm_max_index_all{trial};
    max_inds_all_arenas= max_indices_all{trial};
    rates_all_arenas= peak_rates_all{trial};
    zone_mats_all_arenas=zone_mats_all{trial};
    stability_inds= same_arena_inds{trial}; %arenas with same prefix
    rate_maps_all_arenas= rate_mats_all{trial};
    MD_scores_all_arenas= MD_scores_all{trial};
    gridness_all_arenas=gridness_all{trial};
    autocorrs_all_arenas= autocorrs_all{trial};

    % for rescaling arenas, remove below:
    % check if same context trials is stable enough
    rate_map_1= rate_maps_all_arenas{stability_inds(1)};
    rate_map_2= rate_maps_all_arenas{stability_inds(2)};
    zone_mat_1= zone_mats_all_arenas{stability_inds(1)};
    zone_mat_2= zone_mats_all_arenas{stability_inds(2)};

    [zone_mat_1, zone_mat_2]= ArenaSameSize(zone_mat_1, zone_mat_2);

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

        len=length(rate_maps_all_arenas);

        for a= 1:len-1
            for b=a+1:len

                %if ~isequal([a b], stability_inds) %DONT UNDERSTAND< NOT STABILITY INDICES?
                if true
                    
                rate_mat_1= rate_maps_all_arenas{a};
                rate_mat_2= rate_maps_all_arenas{b};

                max_inds_1= max_inds_all_arenas{a};
                max_inds_2= max_inds_all_arenas{b};

                rates_1=rates_all_arenas{a};
                rates_2=rates_all_arenas{b};

                zone_mat_1= zone_mats_all_arenas{a};
                zone_mat_2= zone_mats_all_arenas{b};

                %                 hyper_ind_1= hyper_inds_all_arenas(a,:);
                %                 hyper_ind_2= hyper_inds_all_arenas(b,:);

                [lenn_1,~] =size(max_inds_1);
                [lenn_2,~]= size(max_inds_2);

                %normalize max inds to arena size
                %                 size_2= size(rate_mat_2);
                %                 norm_max_inds_2=nan(lenn_2,2);
                %                 norm_max_inds_2(:,1)= max_inds_2(:,1)/size_2(1);
                %                 norm_max_inds_2(:,2)= max_inds_2(:,2)/size_2(2);



                %for remapping:
                MD_1= MD_scores_all_arenas(a);
                MD_2= MD_scores_all_arenas(b);

                grid1=gridness_all_arenas(a);
                grid2=gridness_all_arenas(b);

                acorr1= autocorrs_all_arenas{a};
                acorr2= autocorrs_all_arenas{b};

                % if enough number of fields in both trials and all other
                % parameters

                if grid1 >= 0.3 || grid2 >= 0.3

                    if  MD_1 <= 0.25 || MD_2 <= 0.25

                        if lenn_1 >= PF_thresh && lenn_2 >= PF_thresh %number of fields

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

                            spacing1= findPlaceFieldRadius(acorr1) * (10/7);
                            spacing2= findPlaceFieldRadius(acorr2) * (10/7);

                            [phase_shift]= FindPhaseShift(xcorr, [spacing1 spacing2]);

                            xcorr= Cross_Correlation(zone_mat_1,zone_mat_2);
                            R_outer=FindROuter(xcorr);
                            [gridness3] = GridnessRadiusXcorr(xcorr,R_outer);

                            [gc_corr, orig, next]= FindFiringProfileGC(max_inds_1, max_inds_2, rates_1, rates_2);

                            all_corrcoefs(all_count)=gc_corr;
                            phase_shifts(all_count)= phase_shift;
                            xcorr_scores(all_count)=gridness3;
                            all_count=all_count+1;

                            %divide by category

                            if phase_shift < phase_thresh & gridness3 > xcorr_thresh %didnt remap

                                no_remap_corrs(no_remap_count)=gc_corr;
                                no_remap_count=no_remap_count+1;

                            elseif phase_shift >= phase_thresh | gridness3 <= xcorr_thresh

                                remap_corrs(remap_count)=gc_corr;
                                remap_count=remap_count+1;

                            end


                            disp('');
                        end
                        end %of if same context
                    end %of if both arenas pass all criteria
                end % of arena 2
            end %of arena 1
            disp('');
        end %end of if same context phase stable
        %  end % end of if same context xcorr gridness stable
    end % of all cells
end
    cd('C:\Noam\Output\rebekkah\');
    save('remapping vs noremapping corrs no same noam', 'remap_corrs', 'no_remap_corrs',...
        'xcorr_scores','phase_shifts','all_corrcoefs','same_arena_inds')
    cd('C:\Noam\Dropbox\GitTechnion\rebekkah');
% function whatevs
% 
% load('remapping vs noremapping corrs.mat')
% 
% a=remap_corrs;
% b=no_remap_corrs;
% 
% bin_values=-1:0.25:1.1;
% for h=1:length(bin_values)-1
%     x=xcorr_scores;
%     inds=find(x>=bin_values(h) &x<bin_values(h+1));
%     mean_y(h)=nanmean(all_corrcoefs(inds));
%     stderrs(h)=FindStdErr(all_corrcoefs(inds));
%     mean_x(h)=mean(x(x>=bin_values(h) &x<bin_values(h+1)));
% end
% 
% figure;
% 
% subplot(3,3,1)
% PlotMeanErrorbars(a,b);
% box off;
% 
% subplot(3,3,[2 3])
% errorbar(mean_x,mean_y,stderrs);
% box off;
% 
% centers = -0.9:0.1:0.9;
% counts = hist(no_remap_corrs,centers);
% pcts = 100 * counts / sum(counts);
% 
% subplot(3,3,[4 5])
% bar(centers,pcts)
% xlim([-1 1])
% box off;
% 
% centers = -0.9:0.05:0.9;
% counts = hist(remap_corrs,centers);
% pcts = 100 * counts / sum(counts);
% 
% subplot(3,3,[7 8])
% bar(centers,pcts)
% xlim([-1 1])
% box off;
% 
% subplot(3,3,1)
% hold on;
% stderr=FindStdErr(a);
% errorbar(mean(a), stderr, 'o')
% %figure; OverlayingHists(a,b);
% 
% 
% function PlotMeanErrorbars(a,b)
% 
% all_means= [mean(a) mean(b)];
% stderr1=FindStdErr(a);
% stderr2=FindStdErr(b);
% %stderr3=FindStdErr(c);
% all_stderrs= [stderr1 stderr2];
% 
% errorbar(all_means, all_stderrs, 'o')
% 
% function stderr= FindStdErr(x)
% 
% stderr= nanstd(x)/ sqrt((sum(~isnan(x))));