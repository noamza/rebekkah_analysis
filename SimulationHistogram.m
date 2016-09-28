function SimulationHistogram

dbstop if error

parms.dir_load_data = 'N:\users\rebekkah\border effect simulation new';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
    
sums_ns= nan(1,length(file_names));
fanos_med= nan(1,length(file_names));
fanos_ns= nan(1,length(file_names));


count=0;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);
    Cell=dat;
    
   % means_ns(i)= mean(Cell.simulated_distance_ns);
    sums_ns(i)= sum(Cell.simulated_distance_ns <=0.1)/length(Cell.simulated_distance_ns); 
    
   % means(i)= mean(Cell.simulated_distance);
   % sums(i)= sum(Cell.simulated_distance <=0.15); 
    
   % fanos(i)= mean(Cell.simulated_fano_factor);
    fanos_mean(i)= median(Cell.simulated_fano_factor);
    fanos_mean_ns(i)= median(Cell.simulated_fano_factor_ns);
    
   if median(Cell.simulated_fano_factor) > 1
       disp(sprintf('%s', file_name)); 
       
       count=count+1;
   end
   
   
end

count

figure; hist(sums_ns); hold on;
y= 0:50;
x= ones(1,51)* 0.6512; 
plot(x,y,'-');

figure; hist(fanos_mean); hold on;
y= 0:60;
x= ones(1,61)* 2.186; 
plot(x,y,'-');

figure; hist(fanos_mean_ns); hold on;
y= 0:60;
x= ones(1,61)* 2.186; 
plot(x,y,'-');

 disp('');   