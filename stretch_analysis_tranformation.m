parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';
parms.dir_load_data2= 'C:\Users\Dori\Desktop\saved mat with all shapeSeqs';
parms.dir_load_data3= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\cluster matrixs above 0.3';

parms.dir_load_dataVR= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ls to vr\results';
parms.dir_load_dataHR= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ls to hr\results';
parms.dir_load_dataSS= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ls to ss\results';

parms.dir_save_images= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\transformation drawigns one above 0.3';


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

stretch_count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    A=load(file_name);              % opens results for rate mats
    cd(parms.dir_load_data2);
    file_name2=file_name(1:end-10);
    B= load(file_name2);            %opens original file for shapeSeq
    
   cd(parms.dir_load_data3);
    file_name3= sprintf('%s.mat', file_name);
    C= load(file_name3);
    
    cd(parms.dir_load_dataVR);
    
    if exist(file_name, 'file')
        
        VR=load(file_name);
        cd(parms.dir_load_dataHR);
        HR=load(file_name);
        cd(parms.dir_load_dataSS);
        SS=load(file_name);
        
        
        [predicted_SS_stretch]= predictStretchSS(VR.max_stretch, HR.max_stretch);
            
            predicted_minus_actual(i)= predicted_SS_stretch- SS.max_stretch;
            
            random= randsample(1:0.5:5,1); 
            
            random_minus_actual(i)= random - SS.max_stretch;
            
            
    end
end

figure; hist(abs(predicted_minus_actual),3);
figure; hist(abs(random_minus_actual),3);
            