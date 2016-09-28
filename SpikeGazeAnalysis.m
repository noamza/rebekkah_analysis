parms.dir_load_data = 'N:\users\rebekkah\final data smoothed\data sets\results with smoothed gaze mat and disregard when near border [0.05';

dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

for h=1:length(file_names)
    file_name = file_names{h};
    load(file_name);
    
    % find if location of maximum gazing correlates to maximum spiking gaze
    
    [~, gaze_spike_ind] = max(S.gaze_spike_rate(:));
    
    [R,C] = ind2sub(size(S.gaze_spike_rate),gaze_spike_ind);
    gaze_spike_ind(1)=R;
    gaze_spike_ind(2)=C;
    
    [size_x, size_y] = size(S.gaze_rate);
    [value, gaze_ind] = max(S.gaze_rate(:));
    
    [R1,C1] = ind2sub(size(S.gaze_rate),gaze_ind);
    gaze_ind(1)=R1;
    gaze_ind(2)=C1;
    
    if gaze_ind(1) >= 1 && gaze_ind(1) <=4
        gaze_wall = 3; %left wall
    elseif gaze_ind(1) <= size_x && gaze_ind(1) >= size_x-4
        gaze_wall = 4; %right wall
    elseif gaze_ind(2)>= 1 && gaze_ind(2)<= 4
        gaze_wall = 1; %top wall
    elseif gaze_ind(2)<= size_y && gaze_ind(2)>= size_y-4
        gaze_wall =2; %bottom wall
    else
        disp('ERROR ERROR')
    end
    
      if gaze_spike_ind(1) >= 1 && gaze_spike_ind(1) <=4
        gaze_spike_wall = 3; %left wall
    elseif gaze_spike_ind(1) <= size_x && gaze_spike_ind(1) >= size_x-4
        gaze_spike_wall = 4; %right wall
    elseif gaze_spike_ind(2)>= 1 && gaze_spike_ind(2)<= 4
        gaze_spike_wall = 1; %top wall
    elseif gaze_spike_ind(2)<= size_y && gaze_spike_ind(2)>= size_y-4
        gaze_spike_wall =2; %bottom wall
    else
        disp('ERROR ERROR')
      end
    
      if gaze_spike_wall == gaze_wall
          spike_gaze_diff= 0;
      else
          spike_gaze_diff= 2;
      end
      
      spike_gaze_diffs(h)=spike_gaze_diff; 
         
end

      figure; hist(spike_gaze_diffs);    