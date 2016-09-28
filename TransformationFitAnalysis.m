

parms.dir_load_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results above 0.3';
parms.dir_save_data= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\results for transformation from square to smaller rectangle';

parms.dir_save_data2= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\just results square to vertical';

parms.dir_load_data2= 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\cluster matrixs above 0.3'; 

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
    
    %check if first arena is square
    size_arena1= size(A.rate_mats.arena1);
    difference_in_xy= size_arena1(1)-size_arena1(2); %difference between x and y sizes
    
    square=[];
    sm_rectangle=[];
    max_correlation=0;
    max_shift=1;
    max_stretch=1;
    
    if abs(difference_in_xy) < 5       % if arena 1 is a square (plus or minus 5 pixel difference)
        
        square = A.rate_mats.arena1; % will use later only 1st arena shape square
        size_square= size(square);
        
        %use these since rate mat size can vary plus minus 5
        top_limit_x= size_arena1(1)+5;
        low_limit_x= size_arena1(1)-5;
        top_limit_y= size_arena1(2)+5;
        low_limit_y= size_arena1(2)-5;
        
        size_arena2 = size(A.rate_mats.arena2);
        size_arena3 = size(A.rate_mats.arena3);
        size_arena4 = size(A.rate_mats.arena4);
        
        %find the arena with the same x size as square but smaller y axis
        if size_arena2(1)> low_limit_x & size_arena2(1)< top_limit_x & size_arena2(2) < low_limit_y
            sm_rectangle= A.rate_mats.arena2;
        elseif size_arena3(1)> low_limit_x & size_arena3(1)< top_limit_x & size_arena3(2) < low_limit_y
            sm_rectangle= A.rate_mats.arena3;
        elseif size_arena4(1)> low_limit_x & size_arena4(1)< top_limit_x & size_arena4(2) < low_limit_y
            sm_rectangle= A.rate_mats.arena4;
        end
        
        if ~isempty(sm_rectangle)
            
            size_rectangle= size(sm_rectangle);
            
            % whats the best transformation fit for a square arena turning into a smaller rectangle arena [same x, smaller y]?
            
            square_rectangle_diff= size_square(2)-size_rectangle(2);
            stretch_increments= square_rectangle_diff/4;
            
            stretch_count= 0;
            
            correlations= nan([5 5]);
            
            if size_rectangle(1) > size_square(1)
                diff_value= size_rectangle(1)-size_square(1);
                sm_rectangle(size_rectangle(1)-(diff_value-1):size_rectangle(1),:) = [];
            elseif size_rectangle(1) < size_square(1)
                diff_value= size_square(1)-size_rectangle(1);
                square(size_square(1)-(diff_value-1):size_square(1), :) = [];
            end
            
            size_square=size(square);
            size_rectangle=size(sm_rectangle);
            
            
            for stretch_axis= size_rectangle(2):stretch_increments:size_square(2);
                
                factor= stretch_axis/size_rectangle(2);
               % factor= 1+ (1-factor);
                tform = maketform('affine',[factor 0 0; 0 1 0; 0 0 1]);
                %form= [size_square(1), stretch_axis];
                stretch_rect = imtransform(sm_rectangle,tform); % stretches rectangle to square
                
                size_stretch_rect= size(stretch_rect);
                size_diff= size_square(2)- size_rectangle(2);
                
                shift_increments= size_diff/4;
                
                stretch_count=stretch_count+1;
                
                shift_count=0; 
                
                diff_second= size_square(2)- size_stretch_rect(2);
                
%                 if diff_second == 0
%                     
%                     %shift_count= shift_count+1;
%                     shifted_square_portion = square(1:size_stretch_rect(1), 1:size_stretch_rect(2)); 
%                     
%                     correlation= corr2(shifted_square_portion, stretch_rect);
%                     
%                     correlations(stretch_count, shift_count)= correlation;
%                 end  
%                  
                
                for size_skips= 1:shift_increments:diff_second+1;
                    
                    %sizes sometimes off by a bin or two, use larger size
                    %for analysis
%                     if size_square(1)>= size_stretch_rect(1);
%                         use_size = size_stretch_rect(1);
%                     else
%                         use_size= size_square(1);
%                     end
                    
                    shifted_square_portion=nan([size_stretch_rect]);
                    
                    shifted_square_portion(1:size_stretch_rect(1),1:size_stretch_rect(2)) = square(1:size_stretch_rect(1), size_skips:(size_skips+size_stretch_rect(2)-1));
                    
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
            title(sprintf('%0.2f%s %d%s% d%s %0.2f%s',B.cluster_score,'cluster ',max_stretch, 'stretch ', max_shift, 'shift ', max_corr,'max corr '));

            
            subplot(2,3,2);
            imagesc(max_shifted_sq_rm);
            axis equal;
            axis off;
            subplot(2,3,3);
            imagesc(max_stretch_rect_rm);
            axis equal;
            axis off;
             subplot(2,3,5);
            imagesc(square);
            axis equal;
            axis off;
             subplot(2,3,6);
            imagesc(sm_rectangle);
            axis equal;
            axis off;
            title(sprintf('square size=%d%d rectange size=%d%d', size_square, size_rectangle))  
            
            
            
            cd(parms.dir_save_data2)   ;
            title= sprintf('%s',file_name);
            save(title, 'correlations', 'max_shift', 'max_stretch') ;
          %  saveas(fig, sprintf('%s.jpg',file_name))
           % saveas(fig, sprintf('%s.fig',file_name))
            cd(parms.dir_load_data) ;
   
        end
        
        clear correlations; 
    end
    
    clearvars -except file_name i parms file_names 
end