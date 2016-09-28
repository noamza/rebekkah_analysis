parms.dir_load_data = 'N:\users\rebekkah\stability results';
%parms.dir_save_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\cluster criteria 0.45';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

figure;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    peak_rates_begin(peak_rates_begin==0)=nan; 
    peak_rates_begin(peak_rates_end==0)=nan;
    peak_rates_end(peak_rates_end==0)=nan; 
    peak_rates_end(peak_rates_begin==0)=nan; 
    
    scatter(peak_rates_begin, peak_rates_end);
    hold on;
    
    
end