parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

file_names= [{'Results_Cell_r11025_d200505_s01_t6_c2.mat'} ...
    {'Results_Cell_r11340_d221105_s01_t7_c4.mat'} ...
        {'Results_Cell_r13049_d160409_s01,05_t5_c2.mat'}];
    
     % {'Results_Cell_r12138_d300508_sB1_t8_c3.mat'} ...
        
figure; 

count=1;
% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);
    Cell=dat.S;
    
S=[];
simulated_rate_mat=[];
simulated_rate_mat_ns=[];   

strength=1.4;
max_inds= RemoveTooCloseMaxInds(Cell.max_inds, Cell.PF_radius, Cell.rate_mat, strength);

gaussian_mat= createGaussianMat(Cell.rate_mat, Cell.PF_radius, max_inds, mean(Cell.sorted_means));

subplot(3,3,count)
imagesc(Cell.rate_mat); axis square; axis off;
title('Original Rate Map');
subplot(3,3,count+1); 
imagesc(gaussian_mat); axis square; axis off;
title('Generated Equal Rates');

pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;

parms.bin_size=3;
[spk_t, spk_x, spk_y]= Simulate_Spike_Train(pos_mean_x,pos_mean_y,Cell.pos.t,gaussian_mat,parms); 

parms.sigma=1.5;
simulated_rate_map=CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,spk_t,parms);

subplot(3,3,count+2)
imagesc(simulated_rate_map); axis square; axis off;
title('Simulated Spike Trains');

count=count+3;
end

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 17.4])