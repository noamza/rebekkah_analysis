


dbstop if error

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

%dir_name=parms.dir_load_data;
%dir_list = dir(strcat(dir_name,'\*.mat'));
%file_names = {dir_list.name};
count=1;

figure;
n= 2;
m=3;

count=1;

file_names{1}= 'Results_Cell_r11138_d050405_s02_t1_c1.mat';
file_names{2}= 'Results_Cell_r11138_d110405_s01_t5_c1.mat';
file_names{3}= 'Results_Cell_r11138_d130405_s02_t5_c1.mat';
file_names{4}= 'Results_Cell_r11138_d200405_s02_t6_c1.mat';
file_names{5}= 'Results_Cell_r13049_d220309_s01,02,06,07_t6_c3.mat';
file_names{6}= 'Results_Cell_r13049_d220309_s01,02,06,07_t7_c2.mat';


% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;
    
    
    namename= file_name(9:end-4);
    
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
        
        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
     subplot(n,m,count);
            plot(pos_mean_x,pos_mean_y,'k');hold on;
            plot(spk_x,spk_y,'.r');
            axis equal;axis off;
            axis ij
            title(namename, 'fontname', 'calibri', 'Interpreter', 'none', 'fontsize', 10);
            
            subplot(n,m,count+1)
imagesc(Cell.rate_mat);
axis equal; axis off; hold on;
title(sprintf('%0.1f Hz', max(max((Cell.rate_mat)))), 'HorizontalAlignment', 'left', 'fontname', 'calibri', 'fontsize', 10);
plot(Cell.max_index(2), Cell.max_index(1), 'o',  'LineWidth', 4, 'MarkerSize', 20, 'color', 'k'); 

% 'LineWidth', 20,

subplot(n,m,count+2)
imagesc(Cell.zone_mat);
axis equal; axis off;

count=count+3;

end
            