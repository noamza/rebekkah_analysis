parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results 2 above 0.3';
parms.dir_load_data2= 'C:\Users\Dori\Desktop\saved_mat\saved_mat';

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
    
    [size_x, size_y]= size(A.S.rate_mat);
    x_y_diff= abs(size_x-size_y); 
    
    
    
    if isfield(C, 'shapeSeq') & C.shapeSeq{1}== 'ls' | x_y_diff <3    % if arena 1 is a large square
        
        arena1= A.rate_mats.arena1;
        S= get_zones(arena1, A.S);
        arena1_PF_num= S.number_of_PF;  %number of PFs in first arena, large square
        
       
        %% vr arena PF number
        
        index_vr = strfind(C.shapeSeq, 'vr');  %find index of vertical arena
        index_vr=index_vr';
        
        %fill index with nans in empty cells
        ix=cellfun(@isempty,index_vr);
        index_vr(ix)={nan};
        index_vr= find(cell2mat(index_vr)==1);
        
        fldnm = sprintf('arena%d', index_vr);
        rate_mat=[];
        rate_mat = A.rate_mats.(fldnm); %large square 1st arena rate mat
        S=[];
        S=get_zones(rate_mat, A.S);
        vr_PF_num= S.number_of_PF;  % number of PF in vr arena
       
        
        %% hr arena PF number
        
        index_hr = strfind(C.shapeSeq, 'hr');  %find index of vertical arena
        index_hr=index_hr';
        
        %fill index with nans in empty cells
        ix=cellfun(@isempty,index_hr);
        index_hr(ix)={nan};
        index_hr= find(cell2mat(index_hr)==1);
        
        fldnm = sprintf('arena%d', index_hr);
        rate_mat=[];
        rate_mat = A.rate_mats.(fldnm); %large square 1st arena rate mat
        S=[];
        S=get_zones(rate_mat, A.S);
        hr_PF_num= S.number_of_PF;  % number of PF in vr arena
        
         %% hr arena PF number
        
        index_ss = strfind(C.shapeSeq, 'ss');  %find index of vertical arena
        index_ss=index_ss';
        
        %fill index with nans in empty cells
        ix=cellfun(@isempty,index_ss);
        index_ss(ix)={nan};
        index_ss= find(cell2mat(index_ss)==1);
        
        fldnm = sprintf('arena%d', index_ss);
        rate_mat=[];
        rate_mat = A.rate_mats.(fldnm); %large square 1st arena rate mat
        S=[];
        S=get_zones(rate_mat, A.S);
        ss_PF_num= S.number_of_PF;  % number of PF in vr arena
        
       %% ls arena PF number
        rate_mat=[];
        S=[];
       arena5= A.rate_mats.arena5;
        S= get_zones(arena5, A.S);
        arena5_PF_num= S.number_of_PF;  %number of PFs in first arena, large square
      
        %% difference between all types
        
        ls_to_ls(count, 1)= arena1_PF_num;
        ls_to_ls(count, 2)= arena5_PF_num;
        
        ls_to_vr(count, 1)= arena1_PF_num;
        ls_to_vr(count, 2)= vr_PF_num;
        
        ls_to_hr(count, 1)= arena1_PF_num;
        ls_to_hr(count, 2)= hr_PF_num;

        ls_to_ss(count, 1)= arena1_PF_num;
        ls_to_ss(count, 2)= ss_PF_num;
        
        count=count+1;
            
        end
        
end
    
ls_to_all= ls_to_vr;  
ls_to_vr_len= length(ls_to_vr);
ls_to_hr_len= length(ls_to_hr);
ls_to_all(ls_to_vr_len+1:(ls_to_vr_len+ls_to_hr_len), 1)= ls_to_hr(:,1); 
ls_to_all(ls_to_vr_len+1:(ls_to_vr_len+ls_to_hr_len), 2)= ls_to_hr(:,2); 
ls_to_ss_len= length(ls_to_ss);
ls_to_all_len= length(ls_to_all);
ls_to_all(ls_to_all_len+1:(ls_to_all_len+ls_to_ss_len), 1)= ls_to_ss(:,1);
ls_to_all(ls_to_all_len+1:(ls_to_all_len+ls_to_ss_len), 2)= ls_to_ss(:,2);



   % figure; scatter(combo(1,1), combo(1,2), 'o', 'color', 'g')
    
   
   combo= ls_to_vr;
combo(1:15, 3)= ls_to_hr(:,2);
combo(1:15, 4)= ls_to_ss(:,2);


y=0:14;
x=0:14;

    figure; scatter(ls_to_ls(:,1), ls_to_ls(:,2)); title('LS to LS (same arena)'); lsline; hold on; plot(x,y,'-');
    figure; scatter(ls_to_vr(:,1), ls_to_vr(:,2)); title('LS to VR'); lsline; hold on; plot(x,y,'-');
    figure; scatter(ls_to_hr(:,1), ls_to_hr(:,2)); title('LS to HR'); lsline; hold on; plot(x,y,'-');
    figure; scatter(ls_to_ss(:,1), ls_to_ss(:,2)); title('LS to SS'); lsline; hold on; plot(x,y,'-');
    figure; scatter(ls_to_all(:,1), ls_to_all(:,2), 'filled'); title('ALL'); lsline; hold on; plot(x,y,'-');

    

    
    save('place field number changes', 'ls_to_ls', 'ls_to_ss', 'ls_to_vr', 'ls_to_hr');
    