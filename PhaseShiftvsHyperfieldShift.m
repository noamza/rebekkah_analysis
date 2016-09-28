function PhaseShiftvsHyperfieldShift

dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info');

load('remapping data info.mat')

PF_thresh=3; % at least 3 fields
xcorr_thresh= 0;
phase_thresh= 0.2;

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
    zone_mat_1= zone_mats_all_arenas{stability_inds(1)};
    zone_mat_2= zone_mats_all_arenas{stability_inds(2)};
    
    [zone_mat_1, zone_mat_2]= ...
        ArenaSameSize(zone_mat_1, zone_mat_2);
    
    % convert zone mat to 0s and 1s
    zone_mat_1(zone_mat_1~=0)=1;
    zone_mat_2(zone_mat_2~=0)=1;
    
    xcorr= Cross_Correlation(zone_mat_1,zone_mat_2);
    R_outer=FindROuter(xcorr);
    gridness3= GridnessRadiusXcorr(xcorr,R_outer);
    
    if gridness3 > xcorr_thresh % if same context stable enough
        
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
                        (MD_1 <= 0.25 || MD_2 <= 0.25) &&...
                        (grid1>=0.3 || grid2>=0.3)
                    
                    
                  
                    
                    hyperfield_shift=FindDistTwoHypers(hyper_ind_1,hyper_ind_2,norm_max_inds_2);
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
                    
                     if phase_shift >= phase_thresh | gridness3 <= xcorr_thresh
                    
                    hyperfield_shifts(all_count)=hyperfield_shift;
                    all_corrcoefs(all_count)=gc_corr;
                    phase_shifts(all_count)= phase_shift;
                    xcorr_scores(all_count)=gridness3;
                    all_count=all_count+1;
                    
                     end
                end %of if same context
            end %of if both arenas pass all criteria
        end % of arena 2
    end %of arena 1
    disp('');
end %end of if same context phase stable

figure; scatter(phase_shifts, hyperfield_shifts)
lsline

corrcoef(phase_shifts,hyperfield_shifts)

disp('');
