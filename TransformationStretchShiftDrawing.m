
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
        
        ls_size= size(A.rate_mats.arena1);
        
        n=2;
        m=5;
        
        fig=figure;
        subplot(n,m,[1:2 6:7])
        
        x=0;
        y=0;
        w=ls_size(2);
        h=ls_size(1);
        rectangle('Position',[x,y,w,h], 'LineWidth', 3)
        axis ij; axis off; hold on;
        
        
        S= get_zones(A.rate_mats.arena1);
        LS_max_ind= S.max_index;
        
        plot(LS_max_ind(2), LS_max_ind(1), 'x', 'color', 'k', 'MarkerSize', 12, 'LineWidth', 3);

        
        
        subplot(n,m,3)
        imagesc(A.rate_mats.arena1); hold on;
        plot(LS_max_ind(2), LS_max_ind(1), 'x', 'Markersize', 12, 'LineWidth', 3, 'color', 'k') %plots max point in rate mat
        axis equal; axis off; 
        
        %..............................
        
        if isfield(B, 'shapeSeq') & B.shapeSeq{1}== 'ls'     % if arena 1 is a large square
            
            index = strfind(B.shapeSeq, 'vr');  %find index of vertical arena
            index=index';
            
            %fill index with nans in empty cells
            ix=cellfun(@isempty,index);
            index(ix)={nan};
            index= find(cell2mat(index)==1);
            
            fldnm = sprintf('arena%d', index);
            
            vr_size= size(A.rate_mats.(fldnm));
            
            % use stretch to find size
            % use shift to find position
            %plot onto LS
            
            increments= (ls_size(2) - vr_size(2))/4; 
            
             subplot(n,m,[1:2 6:7])
            
            x=0 + (increments * (VR.max_shift-1));
            y=0;
            h=ls_size(1);
            w=vr_size(2) + ((VR.max_stretch-1)* increments);
            rectangle('Position',[x,y,w,h], 'EdgeColor', [1 0 0], 'LineWidth', 3, 'LineStyle', ':');
            axis ij; hold on;
          
            S=[];
            S= get_zones(A.rate_mats.(fldnm));
            VR_max_ind= S.max_index;
            adjust_size= size(VR.max_stretch_rect_rm)./size(A.rate_mats.(fldnm));
            VR_max_ind= VR_max_ind.*adjust_size;
        
            VR_max_ind(2)= VR_max_ind(2) +x;
            VR_max_ind(1)= VR_max_ind(1) +y;
            
             
            plot(VR_max_ind(2), VR_max_ind(1), 'x', 'color', 'r', 'MarkerSize', 12, 'LineWidth', 3);
                   
            
            subplot(n,m,4)
            imagesc(A.rate_mats.(fldnm)); hold on;
            plot(S.max_index(2), S.max_index(1), 'x', 'Markersize', 12, 'LineWidth', 3, 'color', 'k');
            axis equal; axis off;
            title(sprintf('max stretch= %d', VR.max_stretch));
            
            %..........................
            index = strfind(B.shapeSeq, 'hr');  %find index of vertical arena
            index=index';
            
            %fill index with nans in empty cells
            ix=cellfun(@isempty,index);
            index(ix)={nan};
            index= find(cell2mat(index)==1);
            
            fldnm = sprintf('arena%d', index);
            
            hr_size= size(A.rate_mats.(fldnm));
            
            % use stretch to find size
            % use shift to find position
            %plot onto LS
            
            increments= (ls_size(1) - hr_size(1))/4; 
            
             subplot(n,m,[1:2 6:7])
            
            x=0 ;
            y=0 + (increments * (HR.max_shift-1));
            h=hr_size(1) + ((HR.max_stretch-1)* increments);
            w=ls_size(2) ;
            rectangle('Position',[x,y,w,h], 'EdgeColor', [0 1 0], 'LineWidth', 3, 'LineStyle', '--');
            axis ij; hold on;
            
            
            S=[];
            S= get_zones(A.rate_mats.(fldnm));
            HR_max_ind= S.max_index;
             adjust_size= size(HR.max_stretch_rect_rm)./size(A.rate_mats.(fldnm));
            HR_max_ind= HR_max_ind.*adjust_size;
        
            HR_max_ind(2)= HR_max_ind(2) +x;
            HR_max_ind(1)= HR_max_ind(1) +y;
     
            plot(HR_max_ind(2), HR_max_ind(1), 'x', 'color', 'g', 'MarkerSize', 12, 'LineWidth', 3);
            
               subplot(n,m,5)
            imagesc(A.rate_mats.(fldnm)); hold on;
                 plot(S.max_index(2), S.max_index(1), 'x', 'Markersize', 12, 'LineWidth', 3, 'color', 'k');
            axis equal; axis off;
                   title(sprintf('max stretch= %d', HR.max_stretch));
                   
            %..........................
            index = strfind(B.shapeSeq, 'ss');  %find index of vertical arena
            index=index';
            
            %fill index with nans in empty cells
            ix=cellfun(@isempty,index);
            index(ix)={nan};
            index= find(cell2mat(index)==1);
            
            fldnm = sprintf('arena%d', index);
            
            ss_size= size(A.rate_mats.(fldnm));
            
             % use stretch to find size
            % use shift to find position
            %plot onto LS
            
            increments= (ls_size(1) - hr_size(1))/4; 
            
             subplot(n,m,[1:2 6:7])
            
            x=0 + (increments * (SS.max_shift_x-1));
            y=0 + (increments * (SS.max_shift_y-1));
            h=ss_size(1) + ((SS.max_stretch-1)* increments);
            w=ss_size(2) + (increments * (SS.max_stretch-1)) ;
            rectangle('Position',[x,y,w,h], 'EdgeColor', [0 0 1], 'LineWidth', 3, 'LineStyle', '-');
            axis ij; hold on; axis equal;
            
            S=[];
            S= get_zones(A.rate_mats.(fldnm));
            SS_max_ind= S.max_index;
            adjust_size= size(SS.max_stretch_rect_rm)./size(A.rate_mats.(fldnm));
            SS_max_ind= SS_max_ind.*adjust_size;
        
            SS_max_ind(2)= SS_max_ind(2) +x;
            SS_max_ind(1)= SS_max_ind(1) +y;
            
       
            plot(SS_max_ind(2), SS_max_ind(1), 'x', 'color', 'b', 'MarkerSize', 12, 'LineWidth', 3);
            
            
            subplot(n,m,8)
            imagesc(A.rate_mats.(fldnm)); hold on;
                 plot(S.max_index(2), S.max_index(1), 'x', 'Markersize', 12, 'LineWidth', 3, 'color', 'k');
                 axis equal; axis off;
                   title(sprintf('max stretch= %d', SS.max_stretch));
                   
            %..........................

            S= get_zones(A.rate_mats.arena5);
            
            subplot(n,m,9)
            imagesc(A.rate_mats.arena5); hold on;
                 plot(S.max_index(2), S.max_index(1), 'x', 'Markersize', 12, 'LineWidth', 3, 'color', 'k');
                 axis equal; axis off;
            %..........................
            
            
            [predicted_SS_stretch]= predictStretchSS(VR.max_stretch, HR.max_stretch);
            
            predicted_minus_actual(i)= predicted_SS_stretch- SS.max_stretch;
            
            random = randsample(1:5,1); 
            
            subplot(n,m,[1:2 6:7])
            title(sprintf('%0.2f',  C.cluster_score)); 
            
            
            stretch(stretch_count)= VR.max_stretch;
            stretch_count=stretch_count+1;
            stretch(stretch_count)= HR.max_stretch;
            stretch_count=stretch_count+1;
            stretch(stretch_count)= SS.max_stretch;
            stretch_count=stretch_count+1;
            
            cd(parms.dir_save_images)
            saveas(fig, sprintf('%s.jpg', file_name)); 
    
           
            close all
            
        end
    end
end

save('all_stretchs', 'stretch')

figure; hist(stretch);
figure; hist(predicted_minus_actual);