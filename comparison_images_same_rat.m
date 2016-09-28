parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\by rat and arena\17';  
parms.dir_save_pictures = 'N:\users\rebekkah\final data smoothed\comparison images not same date';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('N:\users\rebekkah\final data smoothed\module orientation not same date same rat.mat');


for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};    
    dat = load(file_name); 
    Cell1= dat.S;

    for l= k+1 : length(file_names)
    file_name2 = file_names{l};    
    dat = load(file_name2); 
    Cell2= dat.S;
  
    index= find((module_list(:,1) == Cell1.i) & (module_list(:,2) ==Cell2.i));
    module_diff = module_list(index,3);
    
    fig= figure;
  
    n=2;
    m=3;
    
    fig_pair= subplot(n,m,1);
    fig_pair= plot_ellipse(Cell1.autocorr, Cell1.module.major, Cell1.module.minor ,Cell1.module.phi, Cell1.six_orientation_pts, fig_pair);
    strTitle = sprintf('intersection/score = %0.2f', module_diff);
    title(strTitle);
    
    fig_pair2= subplot(n,m,4);
    fig_pair2= plot_ellipse(Cell2.autocorr, Cell2.module.major, Cell2.module.minor ,Cell2.module.phi, Cell2.six_orientation_pts, fig_pair2);
 
    draw_orientation_line= subplot(n,m,2);
     [draw_orientation_line] = orientation_analysis(Cell1.autocorr, Cell1.rate_mat, Cell1.six_orientation_pts, draw_orientation_line) ;
    
     draw_orientation_line= subplot(n,m,5);
     [draw_orientation_line] = orientation_analysis(Cell2.autocorr, Cell2.rate_mat, Cell2.six_orientation_pts, draw_orientation_line) ;
    
     subplot(n,m,3)
    imagesc(Cell1.zone_mat);
    subplot(n,m,6)
    imagesc(Cell2.zone_mat);
    
    cd(parms.dir_save_pictures);
       %saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
       saveas(fig,sprintf('%d %d.jpg', Cell1.i, Cell2.i'));  
 
       cd(parms.dir_load_data); 
       
       
        close all;
    
    end
end

    