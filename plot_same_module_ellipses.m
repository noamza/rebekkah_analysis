function plot_same_module_ellipses

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\202 cells by rat and arena\2';
parms.dir_save_pictures= 'N:\users\rebekkah\final data smoothed\data sets\202 cells by rat and arena\same module ellipse images';


dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

load('N:\users\rebekkah\final data smoothed\data sets\202 cells by rat and arena\22\module orientation SAME DATE.mat')



for k =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{k};
    dat = load(file_name);
    Cell1= dat.S;
    
    for j=k+1 : length(file_names)
        file_name = file_names{j};
        dat2 = load(file_name);
        Cell2= dat2.S;
        
        for h= 1:length(module_list);
            if Cell1.i == module_list(h,1) & Cell2.i == module_list(h,2) & module_list(h,3)>= 0.7
                fig_pair=figure;
                plot_ellipse_wo_ratemat(Cell1.module.major, Cell1.module.minor ,Cell1.module.phi,Cell1.six_orientation_pts, fig_pair)
                hold on;
                plot_ellipse_wo_ratemat(Cell2.module.major, Cell2.module.minor ,Cell2.module.phi,Cell2.six_orientation_pts, fig_pair)
                
                title(sprintf('angle diff=%.2f, module I/U= %.2f', module_list(h,8), module_list(h,3)));
        
                
                cd(parms.dir_save_pictures);
                %  saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
                saveas(fig_pair,sprintf('%d %d.jpg',Cell1.i,Cell2.i)); %
                cd(parms.dir_load_data);
                
                close all;
            end
        end
        
        
    end
    
end
