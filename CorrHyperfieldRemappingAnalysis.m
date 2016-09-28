
load('zone mats and arena types 5 context remapping.mat')
load('all corrs.mat')

% use cases where at least one arena pair is above 0.49
inds= find(max_corr(max_corr>=0.49)); 

count=1;
for h= 1:length(inds);
    

    trials=inds(h);

    arena_type= arena_types_all{trials};
    zone_mats= zone_mats_all{trials};
    max_indices= norm_max_index_all{trials};
    
    len=5;
    
    for a= 1:len-1 
        for b= a+1:len
            
           % if ~isequal([a b], stab_inds{trials})
            
            zone_mat_1= zone_mats{a};
            zone_mat_2= zone_mats{b};
            
            size_1= size(zone_mat_1);
            size_2= size(zone_mat_2);
            
            new_size= [max(size_1(1),size_2(1)) max(size_1(2),size_2(2))];
            
            max_inds_1= max_indices(a,:);
            max_inds_2= max_indices(b,:);
            
            peak_rates_all_arenas=peak_rates_all{trials};
            
            % stretch so same size if rate mat not exact
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
            
            %find normalized rank of 1st hyperfield to second
            %hyperfield
            % organize rank of second by Dist from first
            
            max_inds_all= max_indices_all{trials};
            max_indices_2= max_inds_all{b};
            [lenn,~]= size(max_indices_2);
            
            max_indices_1= max_inds_all{a};
            peak_rates_1=peak_rates_all_arenas{a};
            ind_1= find(peak_rates_1==max(peak_rates_1));
            max_index_1= max_indices_1(ind_1,:);
            
            field_dists=nan(1,lenn);
            for k=1:lenn
                field_dists(k)= Distance(max_index_1(1),max_index_1(2),...
                    max_indices_2(k,1), max_indices_2(k,2));
            end
            
            % organize by ranked distance
            [~,~,ranks]= unique(field_dists);
            
            [~,p] = sort(ranks,'descend');
            r = 1:length(ranks);
            r(p) = r;
            
            ranks=r/lenn;
            
            peak_rates_2=peak_rates_all_arenas{b};
            max_ind= find(peak_rates_2==max(peak_rates_2));
            rank_num=ranks(max_ind);
            
            arena_type_1= arena_type{a};
            arena_type_2= arena_type{b};
            
            all_dists(count)= Distance(max_inds_1(1), max_inds_1(2), ...
                max_inds_2(1), max_inds_2(2));
            
            rank_dist(count)= rank_num;

            corrs(count)= corr_r;

            count=count+1;
        end
      %  end
    end
end

clear means
clear stderrs

all_dists=rank_dist;

figure; scatter(corrs, all_dists); lsline;


lims= [-0.6:0.3:0.3 0.9];
 
bins=nan(1,length(corrs));
for h=1:length(corrs)
 for k=1:4
    inds=[];
    inds=find(corrs >= lims(k) & corrs < lims(k+1));
    bins(inds)= k;
end
end

for h=1:4
means(h)= mean(all_dists(bins==h)); 
stderrs(h)= std(all_dists(bins==h))/ sqrt(length(all_dists(bins==h))); 
end

figure; errorbar(means, stderrs, 'ro-');

sum(isnan(all_corrcoefs))

disp('');
%figure; scatter(all_corrs, rank_dist);

%save('binned corrs to mean dists UPDATED', 'means', 'stderrs')
