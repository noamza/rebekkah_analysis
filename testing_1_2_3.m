parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';
parms.dir_load_data2= 'C:\Users\Dori\Desktop\saved mat with all shapeSeqs';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    A=load(file_name);
    
    cd(parms.dir_load_data2)
    C=load(file_name(1:end-10));
    
    if ~isfield(C, 'shapeSeq')
        
        disp(sprintf('%s', file_name))
        
    end
    
end

% 214_06-03-10.mat_cell1.mat- ls vr ss hr ls 
% 214_06-03-21.mat_cell1.mat- ls vr hr ss ls
% 214_06-04-21.mat_cell1.mat- ls ls ss vr hr
% 216_06-04-18.mat_cell1.mat- ls hr ss vr ls
% 216_06-05-19.mat_cell3.mat- ls vr hr ss ls
% 216_06-06-22.mat_cell1.mat- ls ss hr vr ls
% 216_06-07-25.mat_cell1.mat- ls vr hr ss ls
% 216_06-08-21.mat_cell1.mat- ls vr ss hr ls
% 217_06-05-18.mat_cell1.mat- ls vr ss hr ls
% 217_06-05-26.mat_cell1.mat- ls vr ss hr ls
% 217_06-07-05.mat_cell1.mat- ls hr vr ss ls
% 217_06-07-31.mat_cell3.mat- ls hr ss vr ls
% 228_06-09-05.mat_cell4.mat- vr hr ls ss vr