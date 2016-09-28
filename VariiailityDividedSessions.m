function VariabilityDividedSessions
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    
    % get all files in folder
    % create file_list
    % create loop
    
    parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\examine first and last 5 mins fanfactor smooth\results';
    
    dir_name= parms.dir_load_data;
    dir_list = dir(strcat(dir_name,'\*.mat'));
    file_names = {dir_list.name};
    
    
    for i=1:length(file_names)
        
        file_name = file_names{i};
        load(file_name);
        
        variances_first_half(i)= var(S.sorted_means_begin)/mean(S.sorted_means_begin);
        
        variances_second_half(i)= var(S.sorted_means_end)/mean(S.sorted_means_end);
        
        variances_total(i) =var(S.sorted_means_total)/mean(S.sorted_means_total);
        
    end
    
    variances_diff= (variances_first_half-variances_second_half);
    
    save('fano factors', 'variances_first_half', 'variances_second_half', 'variances_diff', 'variances_total')
    
    figure; hist(variances_first_half);
    figure; hist(variances_second_half);
    figure; hist(variances_first_half-variances_second_half);