function TestingSecondNormDist

parms.dir_load_data = 'N:\users\rebekkah\bin size 6 nonsmooth\results updated';
parms.dir_save_pictures= 'N:\users\rebekkah\results and info of analysis\final images';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.S;

    size_rate_mat= size(Cell.rate_mat);
    %find distance of second max field location to border
    peak_rates_wo_max= Cell.peak_rates;
    peak_rates_wo_max(peak_rates_wo_max==max(peak_rates_wo_max))= nan;
    second_ind= find(peak_rates_wo_max==max(peak_rates_wo_max));
    second_ind=second_ind(1); %in case same rate take first
    second_max_index= Cell.max_inds(second_ind,:);
    [~,norm_second_dist(i)]= findDistPtToBorder(size_rate_mat(1), size_rate_mat(2), second_max_index);
    
end

sum(norm_second_dist<=0.1)