parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

parms.bin_size=3;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
      
     pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
     pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
    S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=7;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius*(3/7))); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin7(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin7(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
    
     S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=6;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius* (3/6))); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin6(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin6(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
    
     S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=5;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius* (3/5))); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin5(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin5(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
    
     S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=4;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius* (3/4))); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin4(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin4(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
     S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=3;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius)); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin3(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin3(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
    
     S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=2;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius * (3/2))); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin2(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin2(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
    
     S= Cell;
    rate_mat=[];
    zone_mat=[];
        
    parms.bin_size=1;
    
    rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius*3)); 
    
    [size_x,~]= size(rate_mat);
    
    fano_factor_bin1(i)= var(S.sorted_means)/mean(S.sorted_means);
    distances_bin1(i)= S.max_peak_distance/size_x;
    
    %................ no smooth, same bin
    
%     S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%         
%     parms.bin_size=3;
%     
%     rate_mat= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, (S.PF_radius)); 
%     
%     [size_x,~]= size(rate_mat);
%     
%     fano_factor_nonsmooth_sb(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_nonsmooth_sb(i)= S.max_peak_distance/size_x;
%     
%     
%     %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 0.5;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%     
%      [size_x,~]= size(rate_mat);
%     fano_factor_halfsmooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_halfsmooth(i)= S.max_peak_distance/size_x; 
%     
%     
%     %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 1;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%     
%      [size_x,~]= size(rate_mat);
%     fano_factor_onesmooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_onesmooth(i)= S.max_peak_distance/size_x; 
%     
%     %................. full smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.sigma= 1.5;
%      rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%      [size_x,~]= size(rate_mat);
%     fano_factor_smooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_smooth(i)= S.max_peak_distance/size_x; 
%     
%     %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 2;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%     
%      [size_x,~]= size(rate_mat);
%     fano_factor_fullsmooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_fullsmooth(i)= S.max_peak_distance/size_x; 
%     
%         %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 2.5;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%     
%      [size_x,~]= size(rate_mat);
%     fano_factor_2plus_smooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_2plus_smooth(i)= S.max_peak_distance/size_x; 
%     
%         %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 3;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%      [size_x,~]= size(rate_mat);
%     fano_factor_3smooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_3smooth(i)= S.max_peak_distance/size_x; 
%     
%            %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 3.5;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%          [size_x,~]= size(rate_mat);
%     fano_factor_3plus_smooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_3plus_smooth(i)= S.max_peak_distance/size_x; 
%     
%         %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 4;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%          [size_x,~]= size(rate_mat);
%     fano_factor_4smooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_4smooth(i)= S.max_peak_distance/size_x; 
%     
%         %................. half smooth
%      S= Cell;
%     rate_mat=[];
%     zone_mat=[];
%     
%     parms.bin_size=3;
%     parms.sigma= 4.5;
%     rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
%     S= get_zones_nonsmooth(rate_mat, S, S.PF_radius); 
%     
%          [size_x,~]= size(rate_mat);
%     fano_factor_4plus_smooth(i)= var(S.sorted_means)/mean(S.sorted_means);
%     distances_4plus_smooth(i)= S.max_peak_distance/size_x; 
    
end

bins= [1 2 3 4 5 6 7];
fano_factor(1)= sum(fano_factor_bin1>1);
fano_factor(2)= sum(fano_factor_bin2>1);
fano_factor(3)= sum(fano_factor_bin3>1);
fano_factor(4)= sum(fano_factor_bin4>1);
fano_factor(5)= sum(fano_factor_bin5>1);
fano_factor(6)= sum(fano_factor_bin6>1);
fano_factor(7)= sum(fano_factor_bin7>1);

distance(1)= sum(distances_bin1 <0.15);
distance(2)= sum(distances_bin2 <0.15);
distance(3)= sum(distances_bin3< 0.15);
distance(4)= sum(distances_bin4<0.15);
distance(5)= sum(distances_bin5<0.15);
distance(6)= sum(distances_bin6<0.15);
distance(7)= sum(distances_bin7<0.15);

distancemean(1)= mean(distances_bin1);
distancemean(2)= mean(distances_bin2);
distancemean(3)= mean(distances_bin3);
distancemean(4)= mean(distances_bin4);
distancemean(5)= mean(distances_bin5);
distancemean(6)= mean(distances_bin6);
distancemean(7)= mean(distances_bin7);

figure;scatter(bins,fano_factor);
figure; scatter(bins,distance);
figure; scatter(bins,distancemean);

save('var and distance based on bin size nonsmooth', 'fano_factor', 'distance', 'bins') 