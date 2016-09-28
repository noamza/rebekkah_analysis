

parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';
parms.dir_load_data2= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\cluster matrixs above 0.3';
parms.dir_load_data3= 'C:\Users\Dori\Desktop\saved mat with all shapeSeqs';

parms.dir_save_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ls to ss\results';
parms.dir_save_images= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\transformation results\ls to ss\images';


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
    
    if isfield(C, 'shapeSeq') & C.shapeSeq{1}== 'ls'     % if arena 1 is a large square
        
        index = strfind(C.shapeSeq, 'ss');  %find index of vertical arena
        index=index';
        
        %fill index with nans in empty cells
        ix=cellfun(@isempty,index);
        index(ix)={nan};
        
        index= find(cell2mat(index)==1);
        
        fldnm = sprintf('arena%d', index);
        % rectangle = A.rate_mats.(fldnm);
        % square= A.rate_mats.(fldnm);
        % rectangle = A.rate_mats.(fldnm);
        rectangle= A.rate_mats.(fldnm);
        
        % whats the best transformation fit for a small square turning into
        % a vertical smaller rectangle arena?
        
        %   square = A.rate_mats.arena1;
        % rectangle= A.rate_mats.arena1;
        square= A.rate_mats.arena1;
        
        size_square= size(square);
        size_rectangle= size(rectangle);
        
        
        
        
        square_rectangle_diff= size_square(2)-size_rectangle(2);
        stretch_increments= square_rectangle_diff/4;
        
        stretch_count= 0;
        corr_count=1;
        
        
        for stretch_axis= size_rectangle(2):stretch_increments:size_square(2);
            
            factor= stretch_axis/size_rectangle(2);
            % factor= 1+ (1-factor);
            tform = maketform('affine',[factor 0 0; 0 factor 0; 0 0 1]);
            %form= [size_square(1), stretch_axis];
            stretch_rect = imtransform(rectangle,tform); % stretches rectangle to square
            
            stretch_count=stretch_count+1;
            
            size_stretch_rect= size(stretch_rect);
            size_diff_x= size_square(2)- size_rectangle(2);
            size_diff_y= size_square(1)- size_rectangle(1);
            
            shift_increments_x= size_diff_x/4;
            shift_increments_y= size_diff_y/4;
            
            shift_count_x=0;
            
            diff_second_x=[];
            diff_second_x= size_square(2)- size_stretch_rect(2);
            
            diff_second_y=[];
            diff_second_y= size_square(1)- size_stretch_rect(1);
            
            
            
            for size_skips_x= 1:shift_increments_x:diff_second_x+1;
                
                
                shift_count_x= shift_count_x+1;
                
                shift_count_y= 0;
              
                for size_skips_y= 1:shift_increments_y:diff_second_y+1;
                    
                    
                    shifted_square_portion=nan([size_stretch_rect]);
                    
                    shifted_square_portion(1:size_stretch_rect(1),1:size_stretch_rect(2)) = square(size_skips_y:(size_skips_y+size_stretch_rect(1)-1), size_skips_x:(size_skips_x+size_stretch_rect(2)-1));
                    
                    shift_count_y= shift_count_y+1;
                    
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
                    correlations(corr_count)= correlation;

                    corr_count=corr_count+1; 
                    
                    if correlation > max_correlation
                        max_correlation= correlation;
                        max_shifted_sq_rm= shifted_square_portion;
                        max_stretch_rect_rm= stretch_rect;
                        max_stretch= stretch_count;
                        max_shift_x=shift_count_x;
                        max_shift_y=shift_count_y;
                        
                    end
                    
                end
            end
        end
        
        max_corr= max(correlations);
        
        fig=figure;
        subplot(2,3,1);
        %imagesc(correlations);
        title(sprintf('%0.2f%s %d%s %d%s %d%s %0.2f%s',B.cluster_score,'cluster ',max_stretch, 'stretch ', max_shift_x, 'shift_x ', max_shift_y, 'shift_y', max_corr,'max corr '));
        
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
        save(title, 'correlations', 'max_shift_x', 'max_shift_y', 'max_stretch', 'max_stretch_rect_rm') ;
        cd(parms.dir_save_images)
        saveas(fig, sprintf('%s.jpg',file_name))
        saveas(fig, sprintf('%s.fig',file_name))
        cd(parms.dir_load_data) ;
        
        
        close all
    end
    
    clearvars -except file_name i parms file_names file_name2 file_name3
end