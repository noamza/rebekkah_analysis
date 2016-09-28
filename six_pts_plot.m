
parms.dir_load_data1 = 'N:\users\rebekkah\final data smoothed\FINAL ADJUSTED ANGLES';

dir_name= parms.dir_load_data1;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

figure; hold on;

for h =1:length(file_names)
    cd(parms.dir_load_data1);
    file_name = file_names{h};
    load(file_name);
    
    [size_x,size_y]= size(S.autocorr);
    
    normalized_six_pts=S.six_orientation_pts;
    normalized_six_pts(1:6,1)= normalized_six_pts(1:6,1)/size_x;
    normalized_six_pts(1:6,2)= normalized_six_pts(1:6,2)/size_y;
    
    if S.smallest_degree_wall == 1 | S.smallest_degree_wall == 3
        for count=1:6;
            plot(normalized_six_pts(count,1),normalized_six_pts(count,2),'o')
            hold on;
        end
    elseif S.smallest_degree_wall ==2 | S.smallest_degree_wall==4
        for count=1:6;
            plot(normalized_six_pts(count,2),normalized_six_pts(count,1),'o')
            hold on;
        end
    end
    
end
