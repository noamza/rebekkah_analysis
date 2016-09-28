cd('\\192.114.21.198\Dori_Data\data\rebekkah\All rats')

load('zone mats and arena types 5 context remapping.mat')
load('all corrs.mat')

count=1;

%stab_thresh=0.35;

%inds= find(corrs_all>= stab_thresh);

for h= 1:length(corrs_all)   %length(inds);
    
    %trials=inds(h);
    trials=h;
    
    arena_type= arena_types_all{trials};
    zone_mats= zone_mats_all{trials};
    max_indexs= norm_max_index_all{trials};
    max_indices= max_indices_all{trials};
    peak_rates= peak_rates_all{trials};
    rate_mats_all_arenas= rate_mats_all{trials};
    
    len=5;
    
    for a= 1:len-1
        for b= a+1:len
            
            zone_mat_1= zone_mats{a};
            zone_mat_2= zone_mats{b};
            
            size_1= size(zone_mat_1);
            size_2= size(zone_mat_2);
            
            new_size= [max(size_1(1),size_2(1)) max(size_1(2),size_2(2))];
            
            max_inds_1= max_indices{a};
            max_inds_2= max_indices{b};
            
            rates_1= peak_rates{a};
            rates_2= peak_rates{b};
            
            rate_mat_1=rate_mats_all_arenas{a};
            rate_mat_2=rate_mats_all_arenas{b};
            
            if ~isequal(size_1, new_size)
                [zone_mat_1] = ...
                    StretchImage(zone_mat_1, size_1, new_size);
            end
            
            if ~isequal(size_2, new_size)
                [zone_mat_2] = ...
                    StretchImage(zone_mat_2, size_2, new_size);
            end
            
            zone_mat_1(zone_mat_1 ~= 0)= 1;
            zone_mat_2(zone_mat_2 ~= 0)= 1;
            
            corr_r= corr2(zone_mat_1, zone_mat_2);
            
            all_corrs_all_pairs(count)=corr_r;
            
            center_1=size_1/2;
            center_2=size_2/2;
                      
            
            % xcor: acor of the xcorr of the two rate maps
            % 6 peaks of acor of the xcorr
%             xcorr= Cross_Correlation(rate_mat_1,rate_mat_2);
%             auto_max_inds = FindAutoMaxInds(xcorr);
%             
%             xcorr_size=size(xcorr);
%             center= xcorr_size/2;
%             
%             all_dists=nan(1,length(auto_max_inds));
%             for k=1:length(auto_max_inds)
%                all_dists(k)= Distance(center(1),center(2),auto_max_inds(k,1),...
%                    auto_max_inds(k,2));
%             end
%             
%             min_ind=find(all_dists==min(all_dists));
%             min_ind= min_ind(1);
%            % acorr_xcorr= Cross_Correlation(xcorr, xcorr);
%                       
%             phase_shifts(count)= all_dists(min_ind);
%             norm_phase_shifts(count)= all_dists(min_ind)/...
%                 mean(xcorr_size(1),xcorr_size(2));
            %acorr_xcorr_hex(7,:)= [];
            %xcorr_hex = find_six_points(xcorr, 30);
            
            %Pair.xcor=Cross_Correlation(Pair.Cell1.rate_mat,Pair.Cell2.rate_mat);
            
            % computing the new spatial auto correlation of the cross correlation
            %             spatial_xcorr=Pair.xcor;
            %             Pair.acor=Cross_Correlation(spatial_xcorr,spatial_xcorr);
            %             Pair.acor= adjustMapSize(Pair.acor,length(Pair.Cell1.autocorr));
            %
            %             % finding the modul of the auto correlation of the cross correlation
            %             [Pair.acor_of_xcor.module.major,Pair.acor_of_xcor.module.minor,...
            %                 Pair.acor_of_xcor.module.phi,...
            %                 Pair.acor_of_xcor.module.hex_peaks,x0,y0]...
            %                 =Find_Module(Pair.acor,parms);
            %
            %             [ Pair.brillouin_zone,Pair.phase_between,Pair.normDist2,factor]...
            %             =Find_Phase(Pair.acor_of_xcor.module.hex_peaks,...
            %                  Pair.xcor,Pair.acor_of_xcor.module.hex_peaks,...
            %                  Pair.acor_of_xcor.module.phi);
            
            
            % change xcorr size to be odd
            
%             xcorr_len= length(xcorr);
%             
%             acorr_xcorr= adjustMapSize(acorr_xcorr,xcorr_len);
%             
%             try
%                 acorr_xcorr_hex = find_six_points(acorr_xcorr, 30);
%             end
%             
%             if length(acorr_xcorr_hex)==7
%                 
%                 % find acorr_xcorr points, look at closests point in xcorr
%                 % distance of that is phase diff
%                 
%                 [phase_dist]= FindPhaseDiff(acorr_xcorr_hex, xcorr);
%                 
%                 phase_dists(count)=phase_dist/length(acorr_xcorr);
                
                max_index_1= max_indexs(a,:);
                max_index_2= max_indexs(b,:);
                
                hyperfield_dists(count)= Distance(max_index_1(1),max_index_1(2),...
                    max_index_2(1),max_index_2(2));
                
                % organize by ranked distance
      max_inds_all= max_indices_all{trials};
            max_indices_2= max_inds_all{b};
            [lenn,~]= size(max_indices_2);
            
            ind_1= find(rates_1==max(rates_1));
                     
            size_rm=size(zone_mat_2);
            
            norm_max_indices_2= nan(lenn,2);
            norm_max_indices_2(:,1)= max_indices_2(:,1)/size_rm(1);
            norm_max_indices_2(:,2)= max_indices_2(:,2)/size_rm(2);
            
            field_dists=nan(1,lenn);
            for k=1:lenn
                field_dists(k)= Distance(max_index_1(1),max_index_1(2),...
                    norm_max_indices_2(k,1), norm_max_indices_2(k,2));
            end
           
                       
        % subtract minimum distance from actual distance
        % closest field should be at 0

%             % organize by ranked distance
            [~,~,ranks]= unique(field_dists);
            
%             [~,p] = sort(ranks,'descend');
%             r = 1:length(ranks);
%             r(p) = r;
%             
%             ranks=r/lenn;
            
            ranks=ranks/lenn;
            
            max_ind= find(rates_2==max(rates_2));
            rank_num=ranks(max_ind);
                        
             % reverse_rank_dists(count)= rank_num;
                
             rank_dists(count)=rank_num;
             dist_minus_min(count)= hyperfield_dists(count)- min(field_dists);
             
                count=count+1;
                
                % end organize by ranked distance
                
%             else
%                 
%                 phase_dists(count)=nan;
%                 hyperfield_dists(count)=nan;
%                 count=count+1;
%             end
            
            %
            %          [brillouin_zone,phase1,Norm_dist,factor]=...
            %             Find_Phase(acorr_xcorr_hex,xcorr,acorr_xcorr_hex);
            
            h
            count
            
           
            
            disp('');
            
        end
        
    end
end

figure; hist(rank_dists);

%the larger the correlation, the smaller the phase shift (good. as should be.)
% figure;
% scatter(phase_shifts, all_corrs_all_pairs);
% 
% figure;
% scatter(hyperfield_dists, phase_shifts)
% 
% figure;
% scatter(norm_phase_shifts, all_corrs_all_pairs);
% 
% figure;
% scatter(hyperfield_dists, norm_phase_shifts)

save('results remapping', 'all_corrs_all_pairs', 'phase_shifts', ...
    'hyperfield_dists', 'norm_phase_shifts', 'reverse_rank_dists',...
    'dist_minus_min')

 save('remapping results info UPDATED', 'phase_shift_x', 'hyperfield_dists', ...
                'all_corrs_all_pairs', 'phase_shift_y', 'hf_shift_x', ...
                'hf_shift_y');

figure;
scatter(all_corrs_all_pairs,abs(phase_shift_x),'o')
hold on;
scatter(all_corrs_all_pairs, abs(phase_shift_y), 'ro')

figure;
scatter(phase_shift_x, hf_shift_x, 'o')
hold on;
scatter(phase_shift_y, hf_shift_y, 'ro')

            
            

figure;
hist(all_corrs_all_pairs,13)

y=0:120;
x=ones(1,121) * 0.25;

hold on; plot(x,y,'r-')

figure;
scatter(all_corrs_all_pairs, phase_dists);
xlabel('2D correlation')
ylabel('Phase difference')

remap_inds= find(all_corrs_all_pairs < 0.3);
nonremap_inds= find(all_corrs_all_pairs >= 0.3);

figure;
scatter(phase_dists(remap_inds), hyperfield_dists(remap_inds), 'o');
hold on;
scatter(phase_dists(nonremap_inds), hyperfield_dists(nonremap_inds), 'ro');

% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
% hold on;
% hist(diff_rates_corr,30)
% h1 = findobj(gca,'Type','patch');
% set(h1,'facealpha',0.75);
% xlim([0 1])