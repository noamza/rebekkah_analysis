parms.dir_load_data = 'N:\users\rebekkah\bin size 6 nonsmooth\results updated';
parms.dir_load_data2= 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

fano_factorr= nan(1, length(file_names));
fano_factor_ns= nan(1, length(file_names));

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);
    Cell_ns=dat.S;
    cd(parms.dir_load_data2);
    dat2= load(file_name);
    Cell_s= dat2.S;
    
    fano_factorr(i)= var(Cell_s.sorted_means)/mean(Cell_s.sorted_means);
    
    fano_factor_ns(i)= var(Cell_ns.sorted_means)/mean(Cell_ns.sorted_means);
    
end

stderr= std(fano_factorr)/ sqrt(length(fano_factorr))
mean(fano_factorr)
mean(fano_factor_ns)

disp('');