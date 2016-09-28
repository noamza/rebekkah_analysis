function main
parms.dir_load_data=...
    '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\sargolini data regular format'
 %'E:\Pairs of all data\jitter 5 sigma 1 first comes first\Cell not dori';
 parms.dir_save_data=...
'\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\sargolini with histology'
%  parms.dir_save_figures=...
%   'E:\Dori\Desktop\raw data\Cell phase lock and HD histological layers';
% parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
    %midpoint1 (+pi/2),midpoint2(-pi/2)
% parms.num_of_direction_bins=60;
cd('\\192.114.21.198\Dori_Data\data\rebekkah\data sets')
load('layer_cells_list.mat');
cd(parms.dir_load_data);

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

p=[];
HD_ang=[];
HD=[];

phase_lock=[];
phase_HD=[];
index_cell=[];
histo_L=[];
phase_lock_ang=[];
pWrong=[];

 
 
for  i=1:length(file_names)
    
    
      
    i
    file_name = file_names{i};
    dat =load(file_name);

    Cell= dat.S;
    clear dat;
    
    Cell_name=sprintf('Cell_r%s_d%s_s%0.2s_t%d_c%d.mat',Cell.rat,Cell.date,...
       Cell.session,Cell.tetrode,Cell.cell);
   
   
   
%% sargolini Cells   
   ind = strmatch(Cell_name,cell_list);
   
   if ~isempty(ind) && ~isfield(Cell, 'histology')
        Cell.histology_layer=layer_list(ind);
        histo_L(end+1)=layer_list(ind);
        p(end+1)=i;
 %% if dori   
   elseif isfield(Cell, 'histology')
       tmp_layer=Cell.histology;
       layer=update_layer(tmp_layer);
       Cell.histology_layer=layer;
       histo_L(end+1)=layer;
       p(end+1)=i;
   else
       histo_L(end+1)=nan;
       p(end+1)=i;
       pWrong(end+1)=i;
       disp('')
   end
   
%    [Cell.pos.x2,Cell.pos.y2,Cell.spk.x2,Cell.spk.y2,Cell.LFP,Cell.Fs]=...
%        Update_Fields_of_The_Cell(Cell);
   
  %%%%%%%%%%%%%%%%%%% Head Directionality %%%%%%%%%%%%%%%%%%%%%%%%%
%   if ~isempty(Cell.pos.x2)
%   [Cell.rate_ang,Cell.rayleigh,Cell.rayleigh_angle]=Compute_Head_Directionality...
%     (Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2,Cell.spk.x,Cell.spk.y,...
%     Cell.spk.x2,Cell.spk.y2,parms);
%   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase locking  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
%   [Cell.theta,Cell.theta_phase,Cell.phase_locking.pos_phase,Cell.phase_locking.spk_phase]=...
%        Find_Theta_and_Spikes_Phase_and_Pos_Phase(Cell.LFP,Cell.Fs,Cell.spk.t,Cell.pos.t,parms);
%     
% 
%     [Cell.phase_locking.rate_ang,Cell.phase_locking.rayleigh_score,...
%         Cell.phase_locking.rayleigh_angle]=Compute_Phase_Locking...
%              (Cell.phase_locking.spk_phase,Cell.LFP,Cell.Fs,parms);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
 
        S=Cell;
        
    cd(parms.dir_save_data);
    
    save(sprintf('Cell_r%s_d%s_s%s_t%d_c%d.mat',Cell.rat,Cell.date,Cell.session,...
     Cell.tetrode,Cell.cell),'S') ;
   
 cd(parms.dir_load_data);
 
end   
%% populatin  
%  index_cell(end+1)=i;  
%  
%  if isfield(Cell,'rayleigh')
%      HD(end+1)=Cell.rayleigh;
%      HD_ang(end+1)=Cell.rayleigh_angle;
%  else
%      HD(end+1)=nan;
%      HD_ang(end+1)=nan;
%  end
%  
%  
%   if isfield(Cell,'phase_locking')
%      phase_lock(end+1)=Cell.phase_locking.rayleigh_score;
%      phase_lock_ang(end+1)=Cell.phase_locking.rayleigh_angle;
%   else
%       phase_lock(end+1)=nan;
%       phase_lock_ang(end+1)=nan;
%   end
% 
%  if isfield(Cell,'phase_and_HD')
%     phase_HD(end+1)=Cell.phase_and_HD;
%  else
%      phase_HD(end+1)=nan;
%  end
%  
%  
% end
% 
% 
% unique(histo_L)
% 
% figure
% f_size=20;
% plot(phase_lock,HD,'*k','linewidth',5);hold on;
% ylabel('Head Direction','fontsize',f_size);
% xlabel('Phase Locking','fontsize',f_size);
% %title(sprintf('Layer %d',L),'fontsize',f_size);
% ax=0:0.2:1;
% set(gca,'xtick',ax,'fontsize',f_size);
% set(gca,'box','off');
% xlim([0 0.8]);
% ylim([0 0.8]);
% 
% X(:,1)=HD(~isnan(HD));
% X(:,2)=phase_lock(~isnan(HD));
% 
% [COEFF,score] = princomp(X);
% 
% [n,p] = size(X);
% meanX = mean(X,1);
% Xfit = repmat(meanX,n,1) + score(:,1)*COEFF(:,1)';
% plot(Xfit(:,2),Xfit(:,1),'linewidth',3)
% 
% figure;plot(score(:,1),score(:,2),'*k','linewidth',3)
% xlabel('PCA1','fontsize',20);
% ylabel('PCA2','fontsize',20);
% %title('scatter of prefered theta phase Vs prefered HD','fontsize',20);
% ax=get(gca,'xtick');
% set(gca,'xtick',ax,'fontsize',15);



%% sort dori
% HD=HD(157:216);
% phase_lock=phase_lock(157:216);
% histo_L=histo_L(157:216);
%%


% figure
% L=3;
% f_size=20;
% plot(HD(histo_L==L),phase_lock(histo_L==L),...
%     '*r','linewidth',5);
% xlabel('Head Direction','fontsize',f_size);
% ylabel('Phase Locking','fontsize',f_size);
% title(sprintf('Layer %d',L),'fontsize',f_size);
% ax=0:0.2:1;
% set(gca,'xtick',ax,'fontsize',f_size);
% set(gca,'box','off');
% xlim([0 0.8]);
% ylim([0 0.8]);
% 
% 
% figure;hist(phase_lock(histo_L==L));
% xlabel('Phase Locking','fontsize',20);
% ylabel('count','fontsize',20);
% title(sprintf('histogram phase locked layer=%d',L),'fontsize',20)
% 
% 
% figure;hist(HD(histo_L==L));
% xlabel('HD','fontsize',20);
% ylabel('count','fontsize',20);
% title(sprintf('histogram HD layer=%d',L),'fontsize',20)
% 
% 
% 
% 
% figure;hist(HD(histo_L==5));
% %% sargolini Dori Seperation
%  figure;
%  plot(HD(1:156),phase_lock(1:156),'*r','linewidth',5)
%  xlabel('Head Direction','fontsize',20);
%  ylabel('Phase Locking','fontsize',20);
%  title(sprintf('sargolini cells'),'fontsize',20);
%  xlim([0 0.9]);
% ylim([0 0.8]);
% 
%  figure;
%  plot(HD(1:156),phase_lock(1:156),'*r','linewidth',5);hold on;
%  plot(HD(157:216),phase_lock(157:216),'*b','linewidth',5)
%  xlabel('Head Direction','fontsize',20);
%  ylabel('Phase Locking','fontsize',20);
% % title(sprintf('dori cells'),'fontsize',20);
%  legend('Sargolini','Dori');
%  xlim([0 0.9]);
%  ylim([0 0.8]);
%  
%  
 
%  
% 
% disp('')
% 
% 
% 
% function [posx2,posy2,spk_x2,spk_y2,LFP,Fs]=Update_Fields_of_The_Cell(Cell)
% cd('E:\Dori\Desktop\raw data\gridcell_data_sargolini\Data_sargolini');
% 
% %% updating the posx2 and posy2
% file_name=sprintf('%s-%s%s_POS.mat',Cell.rat,Cell.date,Cell.session);
% load(file_name)
% if ~isempty(posx2)
% spk_x2=interp1(post,posx2,Cell.spk.t);
% spk_y2=interp1(post,posy2,Cell.spk.t);
% else
%  posx2=[];
%  posy2=[];
%  spk_x2=[];
%  spk_y2=[];
%  
% end
% 
% %% updating the LFP
% load(sprintf('%s-%s%0.2s_EEG',Cell.rat,Cell.date,Cell.session));
% Fs;
% LFP=EEG;
% 
% if length(Cell.session)>2
%     
%     load(sprintf('%s-%s%s_EEG',Cell.rat,Cell.date,Cell.session(4:5)))
%     LFP=[LFP;EEG];
% end
% 
%     
% disp('')
% 
% 
% 
% function [rate_ang,rayleigh_score,rayleigh_angle]=Compute_Head_Directionality...
%     (pos_x,pos_y,pos_x2,pos_y2,spk_x,spk_y,spk_x2,spk_y2,parms)
% 
% 
% % the head direction of the animal through out the trail
% time_phi = atan2(pos_y2-pos_y,pos_x2-pos_x);
% 
% 
% % the head direction of the animal when spike has ocurred
% count_phi = atan2(spk_y2-spk_y,spk_x2-spk_x);
% 
% n_bins=parms.num_of_direction_bins;
% ang_ax = (-pi+pi/n_bins):pi/(n_bins/2):(pi-pi/n_bins);
% count = hist(count_phi,ang_ax);
% time = hist(time_phi,ang_ax);
% rate_ang=count./time;
% 
% 
% x_val = cos(ang_ax);
% y_val = sin(ang_ax);
% norm_val = nansum(rate_ang);
% x_vec = nansum(cos(ang_ax).*rate_ang);
% y_vec = nansum(sin(ang_ax).*rate_ang);
% vec_len = sqrt(x_vec.^2+y_vec.^2);
% 
% rayleigh_score = vec_len/norm_val;
% rayleigh_angle=rad2deg(atan2(y_vec,x_vec));
% 
% %bar(rad2deg(ang_ax),rate_ang);
% 
% % figure;bar(ang_ax,time); title('time spent in each direction');
% % figure;bar(ang_ax,count);title('number of spikes fired in each direction');
%  
% disp('')
% 
% 
% 
% 
% 

% function layer=update_layer(tmp_layer)
% switch(tmp_layer)
%    
%    
%     case 'MEC'
%        layer=0;
%     case 'L2/3?'
%        layer=0;
%     case 'L2'
%        layer=2;
%     case 'L2?'
%        layer=2;
%     case 'MEC L2'
%        layer=2;
%     case 'MEC II'
%        layer=2;
%     case 'L3'
%        layer=3;
%        case 'L3?'
%        layer=3;
%     case 'EC, L3?'
%        layer=3;
%     case 'EC, L3'
%        layer=3;
%     case 'EC, L3 ?'
%        layer=3;
%     case 'MEC,L3?'
%        layer=3;
%     
%       
% end
% 
% disp('')

% function [theta,theta_phase,pos_phase,spk_phase]=...
%     Find_Theta_and_Spikes_Phase_and_Pos_Phase(LFP,Fs,spk_t,pos_t,parms)
% 
% 
% %% filter theta out of LFP and find theta phase
% theta= Coumpute_theta_band(LFP,Fs);
% hilb_theta=hilbert(theta);
% theta_phase=unwrap(angle(hilb_theta))+parms.beg_cycle;% max point of theta(+0),min ponit (+pi),
%     %midpoint1 (+pi/2),midpoint2(-pi/2)
% 
% bin=1/Fs;
% time=((1:length(theta))*bin)-bin;
% 
% %% find phase of spikes and phase of pos_t
% spk_phase = interp1(time,theta_phase,spk_t,'linear','extrap');
% pos_phase = interp1(time,theta_phase,pos_t,'linear','extrap');
% 
% spk_phase=mod(spk_phase,2*pi);
% pos_phase=mod(pos_phase,2*pi);
% 
% disp('')
% 
% 
% function [theta]= Coumpute_theta_band(LFP,Fs)
% 
%     FFT_LFP= abs(fftshift(fft(LFP)));    
%     FFT_LFP_freq=Fs/2:Fs/(length(LFP)-1):Fs/2;
%       
%     Hd=filt_theta1(Fs); % creat filter
%     
%     theta=filtfilt(Hd.sosMatrix,Hd.ScaleValues,LFP);
%     %FFT_teta=abs(fftshift(fft(gama)));    
%     %FFT_teta_freq=-Fs/2:Fs/(length(gama)-1):Fs/2;
%      
% disp('')
% 
% function Hd=filt_theta1(Fs)
% 
% % Butterworth Bandpass filter designed using FDESIGN.BANDPASS.
% 
% %for details
% % http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html
% 
% % All frequency values are in Hz.
%   % Sampling Frequency
% 
% Fstop1 = 3;           % First Stopband Frequency
% Fpass1 = 4;           % First Passband Frequency
% Fpass2 = 12;           % Second Passband Frequency
% Fstop2 = 15;           % Second Stopband Frequency (needs to be lower than Fs/2 ())
% Astop1 = 30;          % First Stopband Attenuation (dB)
% Apass  = 5;           % Passband Ripple (dB)
% Astop2 = 30;          % Second Stopband Attenuation (dB)
% match  = 'stopband';  % Band to match exactly
% 
% % Construct an FDESIGN object and call its BUTTER method.
% h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
%                       Astop2, Fs);
%                 
% Hd= design(h, 'butter', 'MatchExactly', match);
% disp('')
% 
% function [rate_ang,rayleigh_score,rayleigh_angle]=Compute_Phase_Locking(spk_phase,LFP,Fs,parms)
% 
% %%find theta and theta phase
% theta= Coumpute_theta_band(LFP,Fs);
% hilb_theta=hilbert(theta);
% theta_phase=angle(hilb_theta)+parms.beg_cycle;
% theta_phase=mod(theta_phase,2*pi);
% 
% 
% n_bins=360;
% ang_ax = 0:pi/(n_bins/2):2*pi;
% 
% 
% %find distribution of spike phases
% count = hist(spk_phase,ang_ax);
% 
% % find distribution of phases
% time_bin=1/Fs;
% time = hist(theta_phase,ang_ax)*time_bin;
% 
% rate_ang=count./time;
% 
% 
% x_val = cos(ang_ax);
% y_val = sin(ang_ax);
% norm_val = nansum(rate_ang);
% x_vec = nansum(cos(ang_ax).*rate_ang);
% y_vec = nansum(sin(ang_ax).*rate_ang);
% vec_len = sqrt(x_vec.^2+y_vec.^2);
% 
% rayleigh_score = vec_len/norm_val;
% rayleigh_angle=atan2(y_vec,x_vec);
% %  if rayleigh_angle<0
% %    rayleigh_angle=2*pi+rayleigh_angle;
% %  end
% 
% 
% % 
%   %figure;bar(rad2deg(ang_ax),time); title('time spent in each direction');
%  %figure;bar(rad2deg(ang_ax),count);title('number of spikes fired in each direction');
%  %figure;bar(rad2deg(ang_ax),rate_ang);title('number of spikes fired in each direction');
%   
% disp('')