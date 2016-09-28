function AnalysisSimulations

% creates fake rate maps by creating gaussian functions around max_inds
dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
load('3sets G3MD15PF3 data and results.mat', 'gaussian_mats','PF_radii',...
    'pos_x_all','pos_y_all','pos_t_all', 'zone_mats_all','max_indices')

load('simulated results G3MD15PF3 300x per cell againx3.mat')

parms.bin_size=3;
parms.sigma=1.5;

num=300;

len=length(pos_x_all);

%load('simulated results 367cells 250x per cell.mat')
% all_simulated_dists=nan(len,num);
% all_fano_factors= nan(len,num);
% all_max_over_means= nan(len,num);
% all_stability_corrs=nan(len,num);

for i =359:len

    pos_x=pos_x_all{i};
    pos_y=pos_y_all{i};
    pos_t=pos_t_all{i};
    
    PF_rad=PF_radii(i);
    gaussian_mat=gaussian_mats{i};
    
    zone_mat=zone_mats_all{i};
    max_inds=max_indices{i};
    
    % simulates spike trains using Gaussian map and original trajectory
    [spk_t]= Simulate_Spike_Train(pos_x,pos_y,pos_t,gaussian_mat,parms,num);
    
     simulated_fano= nan(1,num);
     sim_dist= nan(1,num);
     max_over_mean= nan(1,num);
all_corrs=nan(1,num);
     
     
    for h= 1:num
        
        spk_t_use= spk_t(h,:);
        spk_t_use(isnan(spk_t_use))=[];
        
        spk_x=interp1(pos_t,pos_x,spk_t_use);
        spk_y=interp1(pos_t,pos_y,spk_t_use);
        
        % create simulated rate map
        %simulated_rate_map=CreateRateMap(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t_use,parms);
        
        simulated_rate_map= CreateRateMap(pos_x, pos_y, pos_t, spk_x, spk_y, spk_t_use, parms);
     
       % figure; imagesc(simulated_rate_map);
        
        % firing stability of split session:
        [~,all_corrs(h), ~,~,~,~,~,~,~,~,~,~]= ...
    FiringStabilitySplitSessions(pos_x,pos_y,pos_t,spk_t_use,...
    zone_mat, max_inds,parms);
        
        % get all the necessary info
     %   [sim_dist(h), ~]= findHyperDist(simulated_rate_map,PF_rad);
        
        [simulated_fano(h), peak_rates]= findFano(simulated_rate_map,PF_rad);

        % figure; imagesc(simulated_rate_map_ns);
        
      %  sorted_rates= sort(peak_rates);
        
       % sr=sorted_rates;
       % max_over_mean(h)= sr(end)/ mean(sr(1:end-1));
           
        disp('')
        
        % close all

    end
    
    all_fano_factors(i,:)= simulated_fano;
    all_simulated_dists(i,:)=sim_dist;
    all_max_over_means(i,:)= max_over_mean;
    all_stability_corrs(i,:)=all_corrs;
    
    i
  
    save('simulated results G3MD15PF3 300x per cell againx3.mat', ...
        'all_simulated_dists',...)
     'all_fano_factors', 'all_max_over_means','all_stability_corrs')
    
% save('simulated results 367cells 250x per cell.mat', ...
%     'all_fano_factors', ...
%     'all_max_over_means')
        
end
