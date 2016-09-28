function ShuffleMaxPeaks(max_indices,rm_size)

dbstop if error

%cd('N:\users\rebekkah\results and info of analysis');
cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF7')
%load('variability and border distances.mat')
%load('combined three sets results.mat')
load('3sets G3MD25 data and results.mat', 'max_indices_ns','rm_size')

% if max indices under diff name
if ~exist('max_indices', 'var')
max_indices=max_indices_ns;
end

if ~exist('rm_size', 'var')
rm_size=all_rm_sizes;
end

shuffle_times=10000;

% inds5 =[31 32 33 34 35 36 47 48];
% inds=inds5;

% Shuffle

    max_less_pt_05= nan(1,length(shuffle_times));
    max_less_pt_1= nan(1,length(shuffle_times));
    max_less_pt_15= nan(1,length(shuffle_times));
    mean_dist= nan(1,length(shuffle_times));
    second_max_less_pt_1= nan(1,length(shuffle_times));

for k= 1:shuffle_times
    
   % sh_norm_dist= nan(1,length(max_indices(inds5)));
   % sh_second_norm_dist= nan(1,length(max_indices(inds5)));
    
    sh_norm_dist= nan(1,length(max_indices));
    sh_second_norm_dist= nan(1,length(max_indices));
   
    for h= 1:length(max_indices)
        
       % part_h=inds5(h);
       part_h=h; 
       
        % Shuffle max indices
        max_inds= max_indices{part_h};
        order= 1:length(max_inds);
        order= Shuffle(order);
        max_inds= max_inds(order,:);
        
        % Assign a max and second-max
        shuffled_max_ind= max_inds(1,:);
        shuffled_second_max_ind= max_inds(2,:);
        
        [~,sh_norm_dist(h)]= findDistPtToBorder([1 rm_size(part_h,1)], [1 rm_size(part_h,2)], shuffled_max_ind);
        [~,sh_second_norm_dist(h)]= findDistPtToBorder([1 rm_size(h,1)], ...
        [1 rm_size(part_h,2)], shuffled_second_max_ind);
        
    end
    
%     sh_norm_dist(sh_norm_dist<= 0.1)=0;
%     sh_norm_dist(sh_norm_dist> 0.1)=1;
    
    % Count sum that pass criteria for "nearness"
    max_less_pt_05(k)= sum(sh_norm_dist < 0.05)/ length(sh_norm_dist);
    max_less_pt_1(k)= sum(sh_norm_dist < 0.1)/ length(sh_norm_dist);
    max_less_pt_15(k)= sum(sh_norm_dist < 0.15)/ length(sh_norm_dist);
    mean_dist(k)= mean(sh_norm_dist);
    
    second_max_less_pt_1(k)= sum(sh_second_norm_dist < 0.1)/ length(sh_second_norm_dist);
    k
end

save('shuffled max peak distributions', 'max_less_pt_05', ...
    'max_less_pt_1', 'max_less_pt_15', 'mean_dist', 'second_max_less_pt_1'); 