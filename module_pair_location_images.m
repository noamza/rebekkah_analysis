parms.dir_load_data = 'C:\Users\Dori\Desktop\module analysis\17';  
parms.dir_save_pictures = 'C:\Users\Dori\Desktop\bin size 6 nonsmooth\comparison images';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('C:\Users\Dori\Desktop\module analysis\16\module orientation same rat same arena.mat');

modules = 1;

for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};    
    dat = load(file_name); 
    Cell1= dat.S;

    for l= k+1 : length(file_names)
    file_name2 = file_names{l};    
    dat = load(file_name2); 
    Cell2= dat.S;
  
    
    fig= figure;
  
    n=2;
    m=2;
    
    fig_pair= subplot(n,m,1);
    fig_pair= plot_ellipse(Cell1.autocorr, Cell1.module.major, Cell1.module.minor ,Cell1.module.phi, Cell1.six_orientation_pts, fig_pair);
    title('
    
    fig_pair2= subplot(n,m,2);
    fig_pair2= plot_ellipse(Cell2.autocorr, Cell2.module.major, Cell2.module.minor ,Cell2.module.phi, Cell2.six_orientation_pts, fig_pair2);
 
    draw_orientation_line= subplot(n,m,3);
    [draw_orientation_line] = orientation_analysis(Cell1.autocorr, Cell1.rate_mat, Cell1.auto_max_inds, Cell1.six_orientation_pts, draw_orientation_line) ;
    
     draw_orientation_line= subplot(n,m,4);
     [draw_orientation_line] = orientation_analysis(Cell2.autocorr, Cell2.rate_mat, Cell2.auto_max_inds, Cell2.six_orientation_pts, draw_orientation_line) ;
     
    modules = modules+1;
    
    cd(parms.dir_save_pictures);
       %saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
       saveas(fig,sprintf('%d %d.jpg', Cell1.i, Cell2.i'));  
 
       cd(parms.dir_load_data); 
       
       
        close all;
    
    end
end

    