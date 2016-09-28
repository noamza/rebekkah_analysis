parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
%parms.dir_load_data2 = 'C:\Users\Dori\Desktop\saved_mat\saved_mat';

%parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';


cd(parms.dir_load_data);
dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};


 count = 1;
 count2=1;
 
for i=1:length(file_names)
    
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat= load(file_name);
    A=dat.S;
    
   % fanos(i)= var(A.sorted_means)/mean(A.sorted_means);
    
peak_score(i)= A.sorted_means(end)/mean(A.sorted_means(1:end-1)); 

[size_x,~]= size(A.rate_mat);

norm_distances(i)= A.max_peak_distance/size_x; 
%     cd(parms.dir_load_data2)
%     file_name2= file_name(1:end-10);
%     B= load(file_name2);
%     
%     size_r_m= size(A.rate_mats.arena1);
%     
%     arena_r_m(1,count2)= size_r_m(1);
%     arena_r_m(2,count2)= size_r_m(2);
%     
%    
%     
%      if isfield(B, 'shapeSeq') 
%          
%          first_arena(count2)= B.shapeSeq(1); 
%     
%         count2=count2+1; 
%          
% %    max_stretch_vert(count)= A.max_stretch;
% %    
%      else
%          
%          shapeSeq= {'na', 'na', 'na'};
%          
%        first_arena(count2) = shapeSeq(1); 
%    
%        count2=count2+1; 
    
    % end
end

% clear file_names;
% 
% cd(parms.dir_load_data2);
% 
% dir_name= parms.dir_load_data2;
% dir_list = dir(strcat(dir_name,'\*.mat'));
% file_names = {dir_list.name};
% 
% count=1;
% 
% for i=1:length(file_names)
%     
%     file_name = file_names{i};
%     B= load(file_name);
%        
%    max_stretch_horiz(count)= B.max_stretch;
%    
%     count=count+1;
%     
% end
% 
% vert_len= length(max_stretch_vert);
% hor_len= length(max_stretch_horiz);
% max_stretch_all= nan([1 (vert_len+hor_len)]);
% 
% max_stretch_all(1:vert_len)= max_stretch_vert;
% max_stretch_all(vert_len+1:end)= max_stretch_horiz; 
% 
% figure; hist(max_stretch_vert);
% figure; hist(max_stretch_horiz);
% figure; hist(max_stretch_all);
% 
% 
% disp('');
% %pos_t_len


figure; scatter(peak_score, normalized_distances)
xlabel('Max-field firing rate over mean firing rate', 'fontname', 'calibri', 'fontsize', 14)