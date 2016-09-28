dbstop if error

parms.dir_load_data = 'E:\DREADDs again\InjectionExperiment_090215';
parms.dir_save_data = 'N:\users\rebekkah\save dreadds results';

cd(parms.dir_load_data);

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
parms.num_of_direction_bins=120;
parms.bin_size=3;
parms.sigma = 3;
file_names = {dir_list.name};


for i=1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.c;
       
     pos_mean_x=(Cell.S.pos.x1 + Cell.S.pos.x2)/2;
     pos_mean_y=(Cell.S.pos.y1 + Cell.S.pos.y2)/2;

   % spkx=interp1(Cell.S.pos.t,pos_mean_x,Cell.S.spk.t);
  %  spky=interp1(Cell.S.pos.t,pos_mean_y,Cell.S.spk.t);
      
     % Smooth the trajectory
%         Win=hamming(50); % 50 bins (of 0.02sec) are equal to 1 sec
%         Win=Win/sum(Win);
%         x1 = nanconv(Cell.S.pos.x1,Win,'edge','nanout');
%         y1 = nanconv(Cell.S.pos.y1,Win,'edge','nanout');
%         x2 = nanconv(Cell.S.pos.x2,Win,'edge','nanout');
%         y2 = nanconv(Cell.S.pos.y2,Win,'edge','nanout');
% 
%         pos_mean_x = nanmean([x1, x2]');
%         pos_mean_x = pos_mean_x';
%         pos_mean_y = nanmean([y1, y2]');
%         pos_mean_y = pos_mean_y';
  
    rate_mat=CreateRateMap(pos_mean_x,pos_mean_y,Cell.S.pos.t,Cell.S.spk.x,Cell.S.spk.y,Cell.S.spk.t,parms);
    
    figure; imagesc(rate_mat)
    
    peak_rate= max(rate_mat(:));
    mean_rate= nanmean2(rate_mat); 
    
    spkx= Cell.S.spk.x;
    spky= Cell.S.spk.y;
    
    cd(parms.dir_save_data)
    save(file_name, 'peak_rate', 'mean_rate', 'rate_mat', 'pos_mean_x', 'pos_mean_y', 'spkx', 'spky');
    cd(parms.dir_load_data)
    
    figure; plot(pos_mean_x,pos_mean_y,'k'); axis image; hold on;
    plot(Cell.S.spk.x,Cell.S.spk.y,'.r');
    
    clearvars -except parms file_names file_name i
    
end
    
    