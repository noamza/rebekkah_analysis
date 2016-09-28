

parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';
parms.dir_load_data2= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\cluster matrixs above 0.3';
parms.dir_load_data3= 'C:\Users\Dori\Desktop\saved_mat\saved_mat';

parms.dir_save_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ss to hr\results';
parms.dir_save_images= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ss to hr\images';


dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    A=load(file_name);
    cd(parms.dir_load_data2);
    file_name2= sprintf('%s.mat', file_name);
    B= load(file_name2);
    cd(parms.dir_load_data3);
    file_name3=file_name(1:end-10);
    C=load(file_name3);
    
    square=[];
    rectangle=[];
    max_correlation= 0;
    max_shift=1;
    max_strecth=1;
    
    if isfield(C, 'shapeSeq') & C.shapeSeq{1}== 'ss'     % if arena 1 is a large square
        
        index = strfind(C.shapeSeq, 'hr');              % find the horizontal rectangle arena
        index=index';
        
        %fill index with nans in empty cells
        ix=cellfun(@isempty,index);
        index(ix)={nan}; 
      
        index= find(cell2mat(index)==1);
        
        fldnm = sprintf('arena%d', index); 
        rectangle = A.rate_mats.(fldnm);   
        
        % whats the best transformation fit for a small square turning into
        % a vertical smaller rectangle arena?
        
        square = A.rate_mats.arena1; 
        size_square= size(square);
        size_rectangle= size(rectangle);
        
        square_rectangle_diff= size_square(1)-size_rectangle(1);
        stretch_increments= square_rectangle_diff/4;
        
        stretch_count= 0;
        
        correlations= nan([5 5]);
        
        if size_rectangle(2) > size_square(2)
            diff_value= size_rectangle(2)-size_square(2);
            rectangle(:,size_rectangle(2)-(diff_value-1):size_rectangle(2)) = [];
        elseif size_rectangle(2) < size_square(2)
            diff_value= size_square(2)-size_rectangle(2);
            square(:,size_square(2)-(diff_value-1):size_square(2)) = [];
        end
        
        size_square=size(square);
        size_rectangle=size(rectangle);
        
        
        for stretch_axis= size_rectangle(1):stretch_increments:size_square(1);
            
            factor= stretch_axis/size_rectangle(1);
            % factor= 1+ (1-factor);
            tform = maketform('affine',[1 0 0; 0 factor 0; 0 0 1]);
            %form= [size_square(1), stretch_axis];
            stretch_rect = imtransform(rectangle,tform); % stretches rectangle to square
            
            size_stretch_rect= size(stretch_rect);
            size_diff= size_square(1)- size_rectangle(1);
            
            shift_increments= size_diff/4;
            
            stretch_count=stretch_count+1;
            
            shift_count=0;
            
            diff_second= size_square(1)- size_stretch_rect(1);
             
            for size_skips= 1:shift_increments:diff_second+1;
                
                
                shifted_square_portion=nan([size_stretch_rect]);
                
                shifted_square_portion(1:size_stretch_rect(1),1:size_stretch_rect(2)) = square(size_skips:(size_skips+size_stretch_rect(1)-1), 1:size_stretch_rect(2));
                
                shift_count= shift_count+1;
                
                size_portion=size(shifted_square_portion);
                size_stretch_rect=size(stretch_rect);
                
                nan_inds_row=[];
                nan_inds_col=[];
                [nan_inds_row, nan_inds_col]= find(isnan(shifted_square_portion));
                if ~isempty(nan_inds_row)
                    shifted_square_portion(nan_inds_row, nan_inds_col)= 0;
                end
                
                nan_inds_row=[];
                nan_inds_col=[];
                [nan_inds_row, nan_inds_col]= find(isnan(stretch_rect));
                if ~isempty(nan_inds_row)
                    stretch_rect(nan_inds_row, nan_inds_col)= 0;
                end
                
                correlation= corr2(shifted_square_portion, stretch_rect);
                correlations(stretch_count, shift_count)= correlation;
                
                if correlation > max_correlation
                    max_correlation= correlation;
                    max_shifted_sq_rm= shifted_square_portion;
                    max_stretch_rect_rm= stretch_rect;
                    max_stretch= stretch_count;
                    max_shift=shift_count;
                end
                
                
            end
        end
        
        
        fig=figure;
        subplot(2,3,1);
        imagesc(correlations);
        title(sprintf('cluster score=%f, stretch=%d, shift=%d', 'max_corr=%f', B.cluster_score, max_stretch, max_shift, max(correlations)))
        subplot(2,3,2);
        imagesc(max_shifted_sq_rm);
        axis equal;
      %  axis off;
        subplot(2,3,3);
        imagesc(max_stretch_rect_rm);
        axis equal;
       % axis off;
        subplot(2,3,5);
        imagesc(square);
        axis equal;
       % axis off;
        subplot(2,3,6);
        imagesc(rectangle);
        axis equal;
        %axis off;
        title(sprintf('square size=%d%d rectange size=%d%d', size_square, size_rectangle))
        
   
        cd(parms.dir_save_data)   ;
        title= sprintf('%s',file_name);
        save(title, 'correlations', 'max_shift', 'max_stretch', 'max_stretch_rect_rm') ;
        cd(parms.dir_save_images) 
         saveas(fig, sprintf('%s.jpg',file_name))
         saveas(fig, sprintf('%s.fig',file_name))
        cd(parms.dir_load_data) ;
        
    end

clearvars -except file_name i parms file_names file_name2 file_name3
end