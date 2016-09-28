dbstop if error

cd('N:\users\rebekkah\results and info of analysis')
load('rescaling arenas info')

count=1;

for i= 1:length(all_rate_mats)  % open info for every cell
    
    %original arena info
    max_indices= all_max_inds{i};
    max_inds_orig= max_indices{1}; %max inds of first arena
    peak_rates_all_arena= all_peak_rates{i};
    peak_rates_orig= peak_rates_all_arena{1};
    rate_mats_all_arena= all_rate_mats{i};
    rate_mat_orig= rate_mats_all_arena{1};
    zone_mats_all_arena= all_zone_mats{i};
    
    num_of_arenas=length(max_indices);
    
    max_inds_orig_orig= max_inds_orig;
    
    orig_rate_mat= rate_mats_all_arena{1};
    [orig_x orig_y]= size(rate_mats_all_arena{1});
    
    for arena_count= 2:num_of_arenas-1 %for all except last arena
        
        rate_mat_next= rate_mats_all_arena{arena_count};
        
        size_1= length(orig_rate_mat);
        size_2= length(rate_mat_next);
        size_max= max(size_1,size_2); % fi
        
        x_norm_stretch= -1.2;
        max_corr=-inf;
        for A=2.^(-1:0.2:1)
            x_norm_stretch= x_norm_stretch+0.2;
            
            y_norm_stretch= -1.2;
            for C=2.^(-1:0.2:1)  %y stretch
                y_norm_stretch= y_norm_stretch+0.2;
                
                stretch_x= orig_x *A; %size you want rescaled arena to stretch to
                stretch_y= orig_y*C;
                
                [next_x next_y]= size(rate_mat_next); %actual size of rescaled arena
                
                factorx= stretch_x/next_x;
                factory= stretch_y/next_y;
                
                %rescales zone mat to same size as original arena
                tform = maketform('affine',[factory 0 0; 0 factorx 0; 0 0 1]);
                stretch_arena = imtransform(rate_mat_next,tform);
                
                [stretch_x stretch_y]= size(stretch_arena);
                
                diff_x= abs(orig_x - stretch_x);
                diff_y=abs(orig_y - stretch_y);
                
                for B= 1:5:diff_x ;
                    
                    for D= 1:5: diff_y;
                        
                        stretch_arena_final=[];
                        orig_arena_final=[];
                        
                        if orig_x < stretch_x;
                            stretch_arena_final= stretch_arena(B:B+orig_x-1,:);
                            orig_arena_final= orig_rate_mat;
                        elseif orig_x > stretch_x;
                            orig_arena_final= orig_rate_mat(B:B+stretch_x-1,:);
                            stretch_arena_final= stretch_arena;
                        else
                            orig_arena_final=orig_rate_mat;
                            stretch_arena_final=stretch_arena;
                        end
                        
                        if orig_y < stretch_y;
                            stretch_arena_final= stretch_arena_final(:,D:D+orig_y-1);
                            orig_arena_final= orig_arena_final;
                        elseif orig_y > stretch_y;
                            orig_arena_final= orig_arena_final(:,D:D+stretch_y-1);
                            stretch_arena_final= stretch_arena_final;
                        else
                            orig_arena_final=orig_arena_final;
                            stretch_arena_final=stretch_arena_final;
                        end
                        
                        correlation= corr2(orig_arena_final, stretch_arena_final);
                        
                        if correlation > max_corr
                            max_corr(count)= correlation;
                            
                            min_x_stretch(count)=x_norm_stretch;
                            min_y_stretch(count)= y_norm_stretch;
                            
                            A_min(count)=A;
                            C_min(count)=C; 
                        end

                    end
                end
            end
        end
        
        count=count+1;
        
    end %end of arena count
    
    
    
end

count %should be 94
figure; hist(abs(min_x_stretch))
figure; hist(min_y_stretch)
figure; hist(A_min)
figure; hist(C_min)

