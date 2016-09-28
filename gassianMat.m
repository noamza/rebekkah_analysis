% creates fake rate maps by creating gaussian functions around max_inds
% need to determine how to incorperate PF_size into it

dbstop if error

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\3 datasets G3MD25PF3';
parms.dir_save_data= '\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info')
load('All Gaussian Matrixs.mat')
cd(parms.dir_load_data);

%% edited after 11xxx

%for k= 1:20

num=100;

simulated_fano= nan(1,num);
simulated_distance= nan(1,num);
max_over_mean= nan(1,num);
sorted_rates= cell(1,num);

%count=1;

% max_over_mean= nan(1, length(file_names));
% max_over_mean_ns=nan(1,length(file_names));

%figure;
% enumerate on cells

% gaussian_mats=cell(1,length(file_names));

all_fano_factors= cell(1,length(file_names));
% all_simulated_dists=cell(1,length(file_names));
% all_max_over_means= cell(1,length(file_names));
% all_sorted_rates=cell(1,length(file_names));

for i =65:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);
    
    if ~isfield(dat,'S')
    %for Bonnevie data:
        Cell=dat.db.B(1);
        Cell.pos = Cell.pos_data;
        Cell.pos.x=Cell.pos.x1;
        Cell.pos.y=Cell.pos.y1;
        Cell.spk= Cell.spike_data;
        Cell.spk.t= Cell.spk.ts;
    else
        Cell=dat.S;
    end
        
    gaussian_mat=gaussian_mats{i};
    PF_rad= all_PF_radii(i);
    %PF_rad= PF_radius(i);
    
    % strength=1.25;
    % max_inds= RemoveTooCloseMaxInds(Cell.max_inds, Cell.PF_radius, Cell.rate_mat, strength);
    
    % gaussian_mats{i}= createGaussianMat(Cell.rate_mat, Cell.PF_radius, max_inds, mean(Cell.sorted_means));
    
    % end
    
    % save('All Gaussian Matrixs', 'gaussian_mats')
    
    % subplot(3,3,count)
    % imagesc(Cell.rate_mat); axis square; axis off;
    % subplot(3,3,count+1);
    % imagesc(gaussian_mat); axis square; axis off;
    
    % simulates spike trains using Gaussian map and original trajectory
    if ~isempty(Cell.pos.x2)|| ~isfield(Cell.pos,'x2') 
     pos_mean_x=(Cell.pos.x+ Cell.pos.x2)/2;
     pos_mean_y=(Cell.pos.y+ Cell.pos.y2)/2;
    else
       pos_mean_x=(Cell.pos.x); % + Cell.pos.x2)/2;
       pos_mean_y=(Cell.pos.y); % + Cell.pos.y2)/2;
    end
    
    parms.bin_size=3;
    [spk_t]= Simulate_Spike_Train(pos_mean_x,pos_mean_y,Cell.pos.t,gaussian_mat,parms,num);
    
    %%%
    
    for h= 1:num
        
        spk_t_use= spk_t(h,:);
        spk_t_use(isnan(spk_t_use))=[];
        
        spk_x=interp1(Cell.pos.t,pos_mean_x,spk_t_use);
        spk_y=interp1(Cell.pos.t,pos_mean_y,spk_t_use);
        
        % create simulate rate map
        parms.sigma=1.5;
        simulated_rate_map=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,spk_t_use,parms);
        
        % parms.bin_size=6;
        % simulated_rate_map_ns= CreateRateMapNoSmooth(pos_mean_x, pos_mean_y, Cell.pos.t, spk_x, spk_y, spk_t, parms);
        
        % subplot(3,3,count+2)
        % imagesc(simulated_rate_map);
        % axis square; axis off;
        %
        % count=count+3;
        
        %figure; imagesc(simulated_rate_map);
        %figure; imagesc(Cell.rate_mat)
        
        % get all the necessary info
        
%         [simulated_fano(h), simulated_distance(h), peak_rates_all_sm]= ...
%             findFanoAndHyperDist(pos_mean_x,pos_mean_y,Cell.pos.t,...
%             spk_x,spk_y,spk_t_use,simulated_rate_map,PF_rad);
        
        [simulated_fano(h), ~]= findFano(simulated_rate_map,PF_rad);

        % figure; imagesc(simulated_rate_map_ns);
        
       % sorted_rates{h}= sort(peak_rates_all_sm);
        
        %sr=sorted_rates{h};
        %max_over_mean(h)= sr(end)/ mean(sr(1:end-1));
        % max_over_mean_ns(i)= S_ns.sorted_means(end)/ mean(S_ns.sorted_means(1:end-1));
        
        disp('')
        
        close all
        
        
    end
    
    all_fano_factors{i}= simulated_fano;
    all_simulated_dists{i}=simulated_distance;
    all_max_over_means{i}= max_over_mean;
    all_sorted_rates{i}=sorted_rates;
    
    i
    
    cd(parms.dir_save_data);
save('simulated results 367cells 100x per cell cont', ...
    'all_fano_factors', 'all_simulated_dists', ...
    'all_max_over_means', 'all_sorted_rates')
cd(parms.dir_load_data)
        
end

% cd(parms.dir_save_data);
% save('simulated results 60x per cell', ...
%     'all_fano_factors', 'all_simulated_dists', ...
%     'all_max_over_means', 'all_sorted_rates')

% save('simulated results 40x per cell', ...
%     'simulated_fano_factor', 'simulated_distance', ...
%     'max_over_mean', 'sorted_rates')
%
%end