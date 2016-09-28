parms.dir_load_data = 'C:\Users\Dori\Desktop\final data smoothed\results'; 

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for k=1:500;

for i=1:length(file_names)
    file_name = file_names{i};
    load(file_name);
   
    sorted_means= S.sorted_means;
    sorted_means_len=length(S.sorted_means);
   
    remove= randi(sorted_means_len-1,1);
    sorted_means(remove)=[];
    
    coeff_of_var(i)= std(sorted_means)/mean(sorted_means);
    
end

mean_coeff_of_var(k)=mean(coeff_of_var(:));

end

save('average coeff of var', 'mean_coeff_of_var');

disp('');