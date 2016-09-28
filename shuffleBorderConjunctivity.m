dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7');
load('3sets G3MD25 data and results.mat',...
    'peak_rates_ns_all','number_zone_mat','max_indices_ns')

% shuffle_times=10000;
% 
% means_20= nan(1,shuffle_times);
% means_25= nan(1,shuffle_times);
% means_30= nan(1,shuffle_times);
% means_40= nan(1,shuffle_times);
% 
% for k= 1:shuffle_times

    len=length(peak_rates_ns_all);
    
    border_over_total_40= nan(1,len);
    border_over_total_30= nan(1,len);
    border_over_total_25= nan(1,len);
    border_over_total_20= nan(1,len);

for i=1:len
    
    peak_rates= peak_rates_ns_all{i};
    
    max_inds=max_indices_ns{i};
    
    %removes max field:
    max_ind=find(peak_rates==max(peak_rates));
    peak_rates(max_ind)=[];
    max_inds(max_ind,:)= [];
    
    % uncomment below for shuffling:
    %peak_rates= Shuffle(peak_rates);
    sorted_peak_rates= sort(peak_rates);

    [border_over_total_40(i)] = CheckBorderConjunctivity(length(peak_rates),sorted_peak_rates, number_zone_mat{i}, peak_rates, max_inds, 0.4);
    [border_over_total_30(i)] = CheckBorderConjunctivity(length(peak_rates),sorted_peak_rates, number_zone_mat{i}, peak_rates, max_inds,0.3);
    [border_over_total_25(i)] = CheckBorderConjunctivity(length(peak_rates),sorted_peak_rates, number_zone_mat{i}, peak_rates, max_inds,0.25);
    [border_over_total_20(i)] = CheckBorderConjunctivity(length(peak_rates),sorted_peak_rates, number_zone_mat{i}, peak_rates, max_inds,0.2);
    
end


%mean(border_over_total)
%median(border_over_total)
% 
save('border fields in top max firing hyperfield excluded', 'border_over_total_20', 'border_over_total_25',...
    'border_over_total_30', 'border_over_total_40');

% means_20(k)= mean(border_over_total_20);
% means_25(k)= mean(border_over_total_25);
% means_30(k)= mean(border_over_total_30);
% means_40(k)= mean(border_over_total_40);
% 
% k
% 
% end
% 
% save('border conjunctivity shuffled', 'means_40', 'means_30', 'means_20', 'means_25');
% 
% disp('');
