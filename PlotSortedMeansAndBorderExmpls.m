function PlotSortedMeansAndBorderExmpls 

load('variability and border distances UPDATED.mat')

inds= find(norm_dist <0.1 & fano > 2);

inds=inds([1 5 10 11 13 20 28]);

parms.dir_load_data='\\192.114.21.198\Dori_Docs\users\rebekkah\final data smoothed\ROTATED ARENA';
dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count=1;

figure;
n=7;
m=4;

letter= [{'A'} {'B'} {'C'} {'D'} {'E'} {'F'} {'G'} {'H'} {'I'} {'G'}];

%letter= [{'E'} {'F'} {'G'} {'H'} {'I'} {'G'}];

cd(parms.dir_load_data);

%  file_names= [{'Results_Cell_r11138_d121105_s01_t5_c1.mat'} ...
%          {'Results_Cell_r11488_d240506_s01,02,07,08_t2_c5.mat'}...
%             {'Results_Cell_r11938_d181007_s01,02,06,07_t5_c7.mat'}...
%      {'Results_Cell_r11684_d300107_s01,02,07,08_t5_c2.mat'} ...       
%    {'Results_Cell_r11685_d230307_s01,02,07,08_t6_c3.mat'}...
%         {'Results_Cell_r11938_d181007_s01,02,06,07_t8_c1.mat'}...
%          {'Results_Cell_r13049_d100309_s01,03_t6_c1.mat'}];
%          {'Results_Cell_r12138_d300508_sB1_t8_c3.mat'}...
         
          %  ];

      %    {'Results_Cell_r11488_d240506_s01,02,07,08_t2_c5.mat'}...
         %  {'Results_Cell_r13049_d100309_s01,03_t6_c1.mat'}...
         
for f=1:7
    file_name = file_names{inds(f)};
    load(file_name);
    Cell=S;
    
    [~, peak_rates, ~]= getZones(Cell.rate_mat);
    
    % calculate the the rat's head direction (using average of both leds)
    pos_mean_x=(Cell.pos.x + Cell.pos.x2)/2;
    pos_mean_y=(Cell.pos.y + Cell.pos.y2)/2;
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    %spike plot
    subplot(n,m,count)
    plot(pos_mean_x,pos_mean_y,'k');hold on;
    plot(spk_x,spk_y,'.r', 'MarkerSize', 4.5);
    axis image;
    axis off;
    axis ij;
    file= file_name(14:end-4);
    title(sprintf('%s', file), 'Interpreter', 'none');
    text(-0.4,1.05,letter{f},'Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
     
    subplot(n,m, count+1);
    imagesc(Cell.rate_mat); axis image
    axis off; hold on;
    title(sprintf('%0.1f Hz', max(max((Cell.rate_mat)))), 'HorizontalAlignment', 'left');
    plot(Cell.max_index(2), Cell.max_index(1), 'o',  'LineWidth', 3, 'MarkerSize', 20, 'color', 'k');
    
    subplot(n,m,count+2);
    imagesc(Cell.zone_mat); axis image; axis off;
    
    subplot(n,m,count+3);
    plot(1:length(peak_rates), sort(peak_rates), 'ko-', 'MarkerFaceColor', [0.15 0.23 0.37]); 
    box off; grid on; 
    ylabel('firing rate')
    xlim([0 length(peak_rates)+1])
    ylim([0 max(peak_rates)+1])
      
    count=count+4;
    
end


set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4 30.45], 'PaperUnits', 'centimeters', 'PaperSize', [17.4 30.45])

saveas(gcf, 'test.png')

end