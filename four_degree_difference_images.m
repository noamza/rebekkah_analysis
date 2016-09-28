function four_degree_difference_images

parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\202 cells by rat and arena\12';
parms.dir_save_pictures= 'N:\users\rebekkah\final data smoothed\data sets\images 4 angle difference rate and zone mats';


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
            if Cell1.i == module_list(h,1) & Cell2.i == module_list(h,2) & module_list(h,3)>= 0.7...
                    & module_list(h,8) > 3.5 & module_list(h,8) < 5.0
                
                fig_pair=figure;
                subplot(2,2,1)
                imagesc(Cell1.rate_mat);
                subplot(2,2,2)
                imagesc(Cell1.zone_mat);
                subplot(2,2,3)
                imagesc(Cell2.rate_mat)
                subplot(2,2,4)
                imagesc(Cell2.zone_mat)
                
                title(sprintf('angle diff= %f \n angle Cell1 = %f \n angle Cell2= %f', module_list(h,8), Cell1.min_angle, Cell2.min_angle));
                
                cd(parms.dir_save_pictures);
                %  saveas(fig,sprintf('Cell_r%s_d%s_s%s_t%d_c%d.fig',Cell.rat,Cell.date,Cell.session,Cell.tetrode,Cell.cell)); %         % debugger - return
                saveas(fig_pair,sprintf('%d %d.jpg',Cell1.i,Cell2.i)); %
                cd(parms.dir_load_data);
                
                close all;
            end
        end
        
        
    end
    
end
