function HDvsMD

parms.beg_cycle=pi/2;
parms.num_of_direction_bins=120;

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\Bonnevie 2013 original data';

count=1;
[HDs,MDs,count]=find_HD_MD(parms,count);

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\Derdikman Data';

[HDs,MDs,count]=find_HD_MD(parms,count);

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\sargolini with histology';

[HDs,MDs,count]=find_HD_MD(parms,count);

figure; scatter(HDs,MDs);

sr=signrank(HDs,MDs)

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info')
save('HDs and MDs','HDs','MDs')

function [HDs,MDs,count]=find_HD_MD(parms,count)

cd(parms.dir_load_data);

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

% enumerate on cells
for i =1:length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
        
    if ~isfield(dat,'S')
        %for Bonnevie data:
        Cell=dat.db.B(1);
        Cell.pos=Cell.pos_data;
        Cell.pos.x=Cell.pos.x1;
        Cell.pos.y=Cell.pos.y1;
        Cell.spk= Cell.spike_data;
        Cell.spk.t= Cell.spk.ts;
    else
        Cell=dat.S;
    end
    
    if isfield(Cell.pos,'x2')&& ~isempty(Cell.pos.x2)
        
        pos_x=(Cell.pos.x+ Cell.pos.x2)/2;
        pos_y=(Cell.pos.y+ Cell.pos.y2)/2;

       [~,HDs(count),~]=ComputeHeadDirectionalityWithSpeedThreshHold...
    (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2,Cell.spk.t,0,parms);
        
        [~,MDs(count),~]=ComputeMovingDirectionalityWithSpeedThreshHold...
                (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.spk.t,0);
       
            count=count+1;
            
    end
end



