function [intersection,union]=Find_Union(couple)

parms.dir_load_data = 'C:\Users\Dori\Desktop\testing_cells_rebekkah'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

%% try with hex_peaks to see different. change back to 6 orient pts after!!

modules = 1;

for i =1:length(file_names)-1
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell1= dat.S;
    
    for j=i+1 : length(file_names)
        file_name = file_names{j};
        dat = load(file_name);
        Cell2= dat.S;
        
        x10=Cell1.six_orientation_pts(7,1);
        y10=Cell1.six_orientation_pts(7,2);
        beta = Cell1.module.phi ;
        sinbeta = sin(beta);
        cosbeta = cos(beta);
        
        alpha =0:pi/100:2*pi;
        sinalpha = sin(alpha);
        cosalpha = cos(alpha);
        
        x1 = x10 + (Cell1.module.major * cosalpha * cosbeta - Cell1.module.minor * sinalpha * sinbeta);
        y1 = y10 + (Cell1.module.major * cosalpha * sinbeta + Cell1.module.minor * sinalpha * cosbeta);
        
        x20=Cell2.six_orientation_pts(7,1);
        y20=Cell2.six_orientation_pts(7,2);
        beta = Cell2.module.phi ;
        sinbeta = sin(beta);
        cosbeta = cos(beta);
        
        alpha =0:pi/100:2*pi;
        sinalpha = sin(alpha);
        cosalpha = cos(alpha);
        
        x2 = x20 + (Cell2.module.major * cosalpha * cosbeta - Cell2.module.minor * sinalpha * sinbeta);
        y2 = y20 + (Cell2.module.major * cosalpha * sinbeta + Cell2.module.minor * sinalpha * cosbeta);
        
        if isreal(x1) && isreal(x2)&& isreal(y1) && isreal(y2)
            [xa, ya] = polybool('union', x1, y1, x2, y2);
            [xb, yb] = polybool('intersection', x1, y1, x2, y2);
            
            
            intersection = polyarea(xb,yb);
            union=polyarea(xa,ya);
        else
            intersection = nan;
            union=nan;
        end
        
        module_diff = intersection/union;
        
        module_list(modules, 1) = Cell1.i;
        module_list(modules, 2) = Cell2.i;
        module_list(modules, 3) = module_diff;
        
        modules = modules+1;
        
        
        
    end
    
end

save('module list hex peaks updated results', 'module_list');



disp('')