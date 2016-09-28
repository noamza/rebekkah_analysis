function [pos_i, pos_j]= ConvertCoordinates(rate_mat, bin_size, pos_x,pos_y)  
       
        min_pos_x= min(pos_x);
        min_pos_y= min(pos_y);

        [size_x, size_y]=size(rate_mat);
        
            pos_j=round((pos_x-min_pos_x)/ bin_size); 
            pos_i=round((pos_y-min_pos_y)/ bin_size);
             
         pos_i(pos_i<=0)=1;
         pos_j(pos_j<=0)=1;
         pos_i(pos_i>max(size_x))= max(size_x);
         pos_j(pos_j>max(size_y))= max(size_y);