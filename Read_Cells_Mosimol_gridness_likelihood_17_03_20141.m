function main
parms.dir_load_data=...
    'C:\Users\Dori\Desktop\testing top lines'
   % 'E:\Dori\Desktop\raw data\Data_Mosimol\Original data'
    % 'E:\Dori\Desktop\raw data\Pairs of all data\jitter 5 sigma 1 first comes first\ALL CELLS IN DATA'
    %'E:\Dori\Desktop\raw data\Data_Mosimol\Original data'
 
 %parms.dir_save_data=...
%'E:\Dori\Desktop\test';
 %parms.dir_save_figures=...
 % 'E:\Dori\Desktop\raw data\ALL Data 08_04_2014\Cell from smoothed 15 gridness';
parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
    %midpoint1 (+pi/2),midpoint2(-pi/2)
parms.num_of_direction_bins=60;
parms.thresh_pure_conj=0.029;
parms.bin_size=3; % each bin in rate map is 3 cm
parms.time_per_bin=0.02; % pos.t(2)-pos.t(1)
parms.sigma=1.5; % gaussiam smoothing factor of rate map
parms.minRadius=10; 
parms.radiusStep=1;
parms.bin_num=60;
parms.dt=0.02;
parms.pop_size_for_gridness_likelihood=100;
parms.interval_num_hist=30;
parms.likelihood_sig=0.05;

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count=1;
 
  
 
for  i=1:length(file_names)
    
   set(0,'defaulttextinterpreter','none'); 
      
    i
    
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat =load(file_name);
    Cell=dat.S;
%     db= dat.db;
%     
%     Cell.area=db.area;
%     
%     switch Cell.area
%      
%      case 'MEC_II'
%           Cell.histology_layer=2;
%           
%      case 'MEC_III'
%          Cell.histology_layer=3;
%          
%      case 'PaS_D'
%          Cell.histology_layer=7;
%      case  'PaS_S'
%          Cell.histology_layer=8; 
%              
%      case  ''
%          Cell.histology_layer=0; 
%                  
%     end
%      
%   if Cell.histology_layer==2 || Cell.histology_layer==3   
%     
%     Cell.file_name=file_name;
%    Cell.pos.x=db.B(1).pos_data.x1;
%     Cell.pos.y=db.B(1).pos_data.y1;
%     Cell.pos.x2=db.B(1).pos_data.x2;
%     Cell.pos.y2=db.B(1).pos_data.y2;
%     Cell.pos.t=db.B(1).pos_data.t;
      %Cell.pos.offset=db.B(1).pos_data.offset;
%     Cell.spk.t=db.B(1).spike_data.ts;
%     Cell.spk.x=interp1(Cell.pos.t,Cell.pos.x,Cell.spk.t);
%     Cell.spk.y=interp1(Cell.pos.t,Cell.pos.y,Cell.spk.t);
%     Cell.spk.x2=interp1(Cell.pos.t,Cell.pos.x2,Cell.spk.t);
%     Cell.spk.y2=interp1(Cell.pos.t,Cell.pos.y2,Cell.spk.t);
%     Cell.rat=db.rat;
%     Cell.date=db.date;
%     Cell.tetrode=db.tetrode;
%     Cell.cell=db.cell;
%     Cell.directory=db.directory;
%     Cell.file_axona=db.B(1).files;
%     Cell.cut_file=db.B(1).cut_file;
%     Cell.session='B1';
    
    
      %%  plot of all cells

%     if isfield(Cell.spk,'x2')
%     
%     %% Head Directionality
%     [Cell.rate_phi,Cell.rayleigh,Cell.rayleigh_angle]=Compute_Head_Directionality...
%     (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.pos.x2,Cell.pos.y2,Cell.spk.x,...
%     Cell.spk.y,Cell.spk.x2,Cell.spk.y2,parms);
%     pop_HD(i)=Cell.rayleigh;
%     else
%       pop_HD(i)=nan;  
%     end
    
    %% mean pos x y
%     px1(1:length(Cell.pos.x))=Cell.pos.x;
%      px2(1:length(Cell.pos.x))=Cell.pos.x2;
%       py1(1:length(Cell.pos.y))=Cell.pos.y;
%       py2(1:length(Cell.pos.y2))=Cell.pos.y2;
%       
%      pos_x=nanmean([px1;px2]); 
%      pos_y=nanmean([py1;py2]);
%      
%      Win=ones(1,15);
%      Win=Win/sum(Win);
%      
%      pos_y=conv(pos_y,Win,'same');
%      pos_x=conv(pos_x,Win,'same');
%      
%      %%
     
     
     
     
     %spk_x=interp1(Cell.pos.t,pos_x(1:length(Cell.pos.t)),Cell.spk.t);
     %spk_y=interp1(Cell.pos.t,pos_y(1:length(Cell.pos.t)),Cell.spk.t);
     
    Cell=S; 
     rate_mat=Creat_Rate_Map(Cell.pos.x,Cell.pos.y,Cell.pos.t,...
         Cell.spk.x,Cell.spk.y,Cell.spk.t,parms);
     
%    rate_mat=Creat_Rate_Map(pos_x,pos_y,Cell.pos.t,...
 %       spk_x,spk_y,Cell.spk.t,parms);
    
     autocorr=Cross_Correlation...
        (rate_mat,rate_mat);
   
   % figure;
%    imagesc(autocorr)
%     [Cell.module.major,Cell.module.minor,...
%                 Cell.module.phi,...
%                 Cell.hex_peaks,x0,y0]...
%                 =Find_Module(autocorr,parms);
% %             
    Cell.autocorr=Cross_Correlation...
        (Cell.rate_mat,Cell.rate_mat);
 
    figure;
    imagesc(autocorr)
    
    [Cell.module.major,Cell.module.minor,...
                Cell.module.phi,...
                Cell.hex_peaks,x0,y0]...
                =Find_Module(Cell.autocorr,parms);
% %             
     % Calculate the radius that gives the best gridness score
%      R_outer=max(Cell.hex_peaks(:,3))*1.15;

%         R_outer=Find_R_outer(Cell.autocorr,parms);
%         
%         [Cell.gridness1.score,Cell.gridness2.score] = gridnessRadius...
%          (Cell.autocorr,parms,R_outer);
% %     
% %         [Cell.gridness]=Gridness_Likelihood(Cell.pos.x,Cell.pos.y,Cell.pos.t,Cell.spk.t,parms,...
% %             R_outer,Cell.gridness2.score);     
% %           
% %         [Cell.gridness1,Cell.gridness2]=Gridness_Likelihood...
% %             (pos_x,pos_y,Cell.pos.t,Cell.spk.t,parms,...
% %              R_outer,Cell.gridness1.score,Cell.gridness2.score);
%        
%       [Cell.gridness1,Cell.gridness2]=Gridness_Likelihood...
%             (Cell.pos.x,Cell.pos.y,Cell.pos.t,Cell.spk.t,parms,...
%              R_outer,Cell.gridness1.score,Cell.gridness2.score);
%        
%     
%     pop_girdness2(i)=Cell.gridness2.score;
%     
%      pop_girdness2_sig(i)=Cell.gridness2.sig;
%      
%      pop_girdness1(i)=Cell.gridness1.score;
%     
%      pop_girdness1_sig(i)=Cell.gridness1.sig;
    
     
     
    % pop_area{i}=Cell.area;
    % pop_layer=Cell.histology_layer;
     
%      tmpf=cell2mat(db.B(1).files);
%      
%     tmp_rat_files=repmat(db.rat,size(tmpf,1),1);
%     tmp_rat_files(1:size(tmpf,1),7:14)=tmpf;
%     
%     if size(tmpf,1)==2
%     pop_rat_files{count}=tmp_rat_files(1,:);
%     count=count+1;
%      pop_rat_files{count}=tmp_rat_files(2,:);
%      count=count+1;
%     elseif size(tmpf,1)==1
%      pop_rat_files{count}=tmp_rat_files(1,:);
%     count=count+1 ;  
%     end
%      
 fig=Plot_Cell(Cell,parms,i);
    
    
   % S=Cell;
    
    
    
    
        cd(parms.dir_save_data);
        saveas(fig,sprintf('Cell%d_r%s_d%s_t%d_c%d.fig',i,Cell.rat,Cell.date,...
      Cell.tetrode,Cell.cell)) ;
         saveas(fig,sprintf('Cell%d_r%s_d%s_t%d_c%d.jpg',i,Cell.rat,Cell.date,...
      Cell.tetrode,Cell.cell)) ;
%    saveas(fig,sprintf('Cell%d_r%s_d%s_s%s_t%d_c%d.fig',i,Cell.rat,Cell.date,...
%       Cell.session,Cell.tetrode,Cell.cell)) ;
%    saveas(fig,sprintf('Cell%d_r%s_d%s_s%s_t%d_c%d.jpg',i,Cell.rat,Cell.date,...
%       Cell.session,Cell.tetrode,Cell.cell)) ;
%   save(sprintf('Cell%d_r%s_d%s_s%s_t%d_c%d.mat',i,Cell.rat,Cell.date,...
 %     Cell.session,Cell.tetrode,Cell.cell),'S');
  
 
    
    
    
    clear dat S Cell;
    close all;
 %end
end

figure;
plot(pop_girdness2,pop_HD,'*k')

sum(pop_girdness2>0.35)

figure;
hist(pop_girdness2,17)


%k=unique(pop_rat_files);

for i=1:length(k)
   tmp=k{i};
   
   rat(i,:)=tmp(1:5);
   
   
   tmp_f=num2str(tmp(7:14));
   tmp_f(9:12)='.EEG';
   file(i,:)=tmp_f;
   
   clear tmp_f
   
end


disp('')


function [gridness1,gridness2]=Gridness_Likelihood(pos_x,pos_y,pos_t,spk_t,...
        parms,R_outer,score1,score2)
    
gridness2.score=score2;
gridness1.score=score1;
    n=parms.pop_size_for_gridness_likelihood;
for i=1:n
    
    if mod(i,10) == 0
      fprintf('%d ',i);
    end
    
    
    
    new_spkt=Shift_Spike_Times(spk_t,pos_t);
    
    % adujst the spike posiotion so they will fit the rat path
    spkx=interp1(pos_t,pos_x(1:length(pos_t)),new_spkt);
    spky=interp1(pos_t,pos_y(1:length(pos_t)),new_spkt);
    
    
    % generate new rate map
    map = Creat_Rate_Map(pos_x,pos_y,pos_t,...
        spkx,spky,new_spkt,parms);
    
    
    
      
   
    autocorr1=Cross_Correlation...
        (map,map);
 
    
       
    
    %%
    [gridness1.pop(i),gridness2.pop(i)] = gridnessRadius...
        (autocorr1,parms,R_outer);
     
        
    
   
end



gridness2.sig=1-sum(gridness2.score>...
    gridness2.pop)/parms.pop_size_for_gridness_likelihood;

gridness1.sig=1-sum(gridness1.score>...
    gridness1.pop)/parms.pop_size_for_gridness_likelihood;

%figure;
%imagesc(Cell.rate_mat)
disp('')

function new_spkt=Shift_Spike_Times(spkt,post)


rnd_shift=max(post)*rand(1);
tmp=spkt+rnd_shift;
ind=find(tmp>max(post));
if length(ind)>0
    new_spkt=tmp(ind)-max(post);
    new_spkt(length(ind)+1:length(ind)+ind(1)-1)=tmp(1:ind(1)-1);
else
    new_spkt=tmp;
end
disp('')

function [fig_Pair]=Plot_Spatial_AutoCorrelation(autocorr,major,minor,phi,hex_peaks,fig_Pair,m,n,l)
    subplot(m,n,l); 
    lw=0.5;
    h=imagesc(autocorr);axis('xy'); axis equal;axis off;hold on;

%center point
x0=hex_peaks(7,1); 
y0=hex_peaks(7,2);
 beta =phi ;
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha =0:pi/100:2*pi;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    x1 = x0 + (major * cosalpha * cosbeta - minor * sinalpha * sinbeta);
    y1 = y0 + (major * cosalpha * sinbeta + minor * sinalpha * cosbeta);

 plot(x1,y1,'k','LineWidth',lw);hold on;
 plot(hex_peaks(:,1),hex_peaks(:,2),'ok','LineWidth',lw);hold on;
 
 xMajor1 = x0 + major * cos(phi);
xMajor2 = x0 - major * cos(phi);
yMajor1 = y0 + major * sin(phi);
yMajor2 = y0 - major * sin(phi);
 p1=xMajor1:(xMajor2-xMajor1)/10:xMajor2;
  p2=yMajor1:(yMajor2-yMajor1)/10:yMajor2;
  if ~isempty(p1) & ~isempty(p2)
plot(p1,p2,'LineWidth',lw); hold on;
  end
 xMinor1 = x0 + minor * cos(phi+pi/2);
xMinor2 = x0 - minor * cos(phi+pi/2);
yMinor1 = y0 + minor * sin(phi+pi/2);
yMinor2 = y0 - minor * sin(phi+pi/2);

p11=xMinor1:(xMinor2-xMinor1)/10:xMinor2;
p21=yMinor1:(yMinor2-yMinor1)/10:yMinor2;
if ~isempty(p11) & ~isempty(p21)
plot(p11,p21,'k','LineWidth',lw); hold on;
end
disp('')

function fig=Plot_Cell(Cell,parms,i)
     
    fig=figure;
    
    m=2;
    n=3;
    count=1;
    
    
    subplot(m,n,count);
    plot(Cell.pos.x,Cell.pos.y,'k',Cell.spk.x,Cell.spk.y,'.r');
    axis equal;axis off;
   
    count=count+1;
    
    subplot(m,n,count);
    imagesc(Cell.rate_mat);axis('xy');axis equal;axis off;
    
    count=count+1;
    
    subplot(m,n,count);
    [fig]=Plot_Spatial_AutoCorrelation...
        (Cell.autocorr,Cell.module.major,Cell.module.minor,...
        Cell.module.phi,Cell.hex_peaks,fig,m,n,count);
%   imagesc(Cell.autocorr);axis('xy');axis equal;axis off;
    count=count+1;
    
    subplot(m,n,count);
    bar(rad2deg(-pi:2*pi/(length(Cell.rate_phi)-1):pi),Cell.rate_phi);
    set(gca,'xtick',-180:120:180);
    title(sprintf('HD Rayleigh=%0.3f\nHD pref=%0.3f',...
    Cell.rayleigh,Cell.rayleigh_angle))
    
    count=count+1;
    subplot(m,n,count);
%     text(0.2,0.8,sprintf('Cell=%d\nlayer=%d\nrat=%s\ndate=%s\ntetrode=%d\ncell=%d\ngridness1=%0.3f\ngridness1 sig=%0.3f\ngridness2=%0.3f\ngridness2 sig=%0.3f'...
%         ,i,Cell.histology_layer,Cell.rat...
%         ,Cell.date,Cell.tetrode,Cell.cell,Cell.gridness1.score,Cell.gridness1.sig,...
%         Cell.gridness2.score,Cell.gridness2.sig));
   % text(0.2,0.8,sprintf('Cell=%d\nrat=%s\ndate=%s\ntetrode=%d\ncell=%d\ngridness1=%0.3f\ngridness1 sig=%0.3f\ngridness2=%0.3f\ngridness2 sig=%0.3f'...
    %    ,i,Cell.rat...
     %   ,Cell.date,Cell.tetrode,Cell.cell,Cell.gridness1.score,Cell.gridness1.sig,...
      %  Cell.gridness2.score,Cell.gridness2.sig));
        axis off
        axis off
    count=count+1;
    
    subplot(m,n,count);
    text(0.2,0.9,sprintf('Module\nmajor=%0.3f\nminor=%0.3f\nphi=%0.3f',...
        Cell.module.major,Cell.module.minor,...
        rad2deg(Cell.module.phi)));
        axis off
    count=count+1;
    
    
        
disp('')
        




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated the position angles according to the angle tAngle [radians]
function [newX,newY] = rotatePath(x,y,tAngle)

newX = x * cos(tAngle) - y * sin(tAngle); % Tilted x-coord
newY = x * sin(tAngle) + y * cos(tAngle); % Tilted y-coord


function [gridness1,gridness2] = gridnessRadius...
        (org_Acor,parms,R_outer)

    
   %% compute distance from center
[val,center_x]=max(max(org_Acor));
[val,center_y]=max(max(org_Acor'));
[Y,X] = ndgrid(1:1:size(org_Acor,1), 1:1:size(org_Acor,2));
dist_from_center=sqrt((Y-center_y).^2+(X-center_x).^2);


%% making sure that outer radius of the ring (R_outer) is not bigger than the distance matrix
R_outer=min([min(dist_from_center(1,:)),min(dist_from_center(:,1)),...
    min(dist_from_center(size(dist_from_center,1),:))...
min(dist_from_center(:,size(dist_from_center,2))),R_outer]);

%% compute inner radius of the anulus (ring)
R_inner=ceil(min(dist_from_center(org_Acor<0.1)));

    
%extract the original anulus (ring) from Acor
org_Ring=org_Acor(dist_from_center<=R_outer & dist_from_center>=R_inner);


%% plot original ring
%  figure;
%  tmp=org_Acor;
% tmp(dist_from_center>R_outer | dist_from_center<R_inner)=nan
%  imagesc(tmp)
 

% make sure that after rotation and interpulation the center will remain the maximum point.
    org_Acor(center_y,center_x)=10;

    for jj = 2:6
         
       %% rotate the auto-correlation 
       rot_Acor=imrotate(org_Acor,(jj-1)*30,'bicubic');
                              
         %% compute distance from new center
        [val,tmp_center_x]=max(max(rot_Acor));
        [val,tmp_center_y]=max(max(rot_Acor'));
        [Y,X] = ndgrid(1:1:size(rot_Acor,1), 1:1:size(rot_Acor,2));
        tmp_dist_from_center=sqrt((Y-tmp_center_y).^2+(X-tmp_center_x).^2);
        
        % extract the anulus(ring)
        rot_Ring=rot_Acor(tmp_dist_from_center<=R_outer & tmp_dist_from_center>=R_inner);
        
%         figure;
%          tmp_rot1=rot_Acor
%          tmp_rot1(tmp_dist_from_center>R_outer | tmp_dist_from_center<R_inner)=nan
%           imagesc(tmp_rot1)
%          
        if length(rot_Ring)~=length(org_Ring) 
            gridness2=nan;
            return 
        end
        
        %% compute pearson correlation between rotate Acor and original Acor
        corrValues(jj) = PointCorr(org_Ring,rot_Ring);
        
       clear rot_Ring tmp_center_x tmp_center_y tmp_dist_from_center Y X
        
    end
    
     

%    % min of higher correlation at 60 and 120 degree rotation
     min_rot_60_120 = min([corrValues(3),corrValues(5)]);
    % max of lower correlations 30,90,150 rotation
    max_rot_30_90_150 = max([corrValues(2),corrValues(4),corrValues(6)]);
    
    %% calculate gridness min(60,120)-max(30,90,150)
    gridness2 = min_rot_60_120 - max_rot_30_90_150;

%% different way to calculate gridness
    gridness1=mean(([corrValues(3),corrValues(5)]))-...
         mean([corrValues(2),corrValues(4),corrValues(6)]);
disp('')

function [rate_ang,rayleigh_score,rayleigh_angle]=Compute_Head_Directionality...
    (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_x,spk_y,spk_x2,spk_y2,parms)

dt=pos_t(2)-pos_t(1);
% the head direction of the animal through out the trail
time_phi = atan2(pos_y2-pos_y,pos_x2-pos_x);


% the head direction of the animal when spike has ocurred
count_phi = atan2(spk_y2-spk_y,spk_x2-spk_x);

step_len=deg2rad(3);%parms.num_of_direction_bins;
ang_ax = -pi:step_len:pi;
count = hist(count_phi,ang_ax);
time = hist(time_phi,ang_ax)*dt;
tmp_rate_ang=count./time;

Win=hamming(10);
Win=Win/sum(Win);

tmp_rate_ang=cconv(tmp_rate_ang,Win');
rate_ang=tmp_rate_ang((length(Win)/2):length(tmp_rate_ang)-(length(Win)/2));


x_val = cos(ang_ax);
y_val = sin(ang_ax);
norm_val = nansum(rate_ang);
x_vec = nansum(cos(ang_ax).*rate_ang);
y_vec = nansum(sin(ang_ax).*rate_ang);
vec_len = sqrt(x_vec.^2+y_vec.^2);

rayleigh_score = vec_len/norm_val;
rayleigh_angle=rad2deg(atan2(y_vec,x_vec));

%bar(rad2deg(ang_ax),rate_ang);

%  figure;bar(ang_ax,time); title('time spent in each direction');
%  figure;bar(ang_ax,count);title('number of spikes fired in each direction');
%  figure;bar(ang_ax,rate_ang);title('number of spikes fired in each direction');
%  
disp('')

function rate_mat=Creat_Rate_Map(posx,posy,post,spkx,spky,spkt,parms)
max_x = (max(posx)); 
max_y = (max(posy));
min_x = min((posx));
min_y = min((posy));

% divid the environment into spatial bins 
axis_x = min_x:parms.bin_size:max_x;
axis_y = min_y:parms.bin_size:max_y;
dt=post(2)-post(1);

time_mat = zeros(length(axis_y),length(axis_x));
spike_mat = zeros(length(axis_y),length(axis_x));
rate_mat = zeros(length(axis_y),length(axis_x));

%create time mat (compute how much time the rat spends in each bin)
% find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to
% 
for i = 1:length(post)
    if ~isnan(posx(i)) && ~isnan(posy(i))
        [min_val,x_ind] =  min(abs(posx(i)-axis_x));
        [min_val,y_ind] =  min(abs(posy(i)-axis_y));
        time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind)+dt;
        
    end
end
%create conut mat( count the num of spikes in each bin)
for i = 1:length(spkt)
   
        [min_val,x_ind] =  min(abs(spkx(i)-axis_x));
        [min_val,y_ind] =  min(abs(spky(i)-axis_y));
        spike_mat(y_ind,x_ind)= spike_mat(y_ind,x_ind)+1;
end

%
%spike_mat=Smooth_Rate_Mat(spike_mat,parms);

%time_mat=Smooth_Rate_Mat(time_mat,parms);
% create rate mat
 rate_mat=spike_mat./time_mat;
 rate_mat(rate_mat==inf)=NaN;
 rate_mat=Smooth_Rate_Mat(rate_mat,parms);
 
%   figure;imagesc(rate_mat);axis('xy')
%   
%   figure;
%   plot(posx,posy,'k',spkx,spky,'.r')
disp('')

function rate_mat=Smooth_Rate_Mat(rate_mat,parms) 
 %create window
 sigma=parms.sigma;
 size_h=[13 13]; %7*[floor(parms.sigma),floor(parms.sigma)];
  h=fspecial('gaussian',size_h,sigma);
  
  %figure;imagesc(h)
 %smooth the rate mat
rate_mat=nanconv2(rate_mat,h);
%figure;imagesc(rate_mat);title('rate map');axis equal
function out_mat=Cross_Correlation(mat1,mat2)

 [ma,na] = size(mat1);
 [mb,nb] = size(mat2);
 mc = max([ma+mb-1,ma,mb]);
 nc = max([na+nb-1,na,nb]);
 m=min(mc,nc);
 out_mat = nan(m,m);
 

i_size = size(mat2,1); j_size = size(mat2,2);
[work_mat,npad_i,npad_j] = pad_edges(mat1,mat2,1);

   for i = 1:size(out_mat,1)
        for j = 1:size(out_mat,2)
                
        % for each i and j, choose the correct sub-mat (size of mat 2) to
        % multiply with mat2
        
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
                                                          npad_j+j-floor(j_size):npad_j+j-1  ); 
        nan_sub_mat=sub_mat .* mat2;                                             
        notnan_inds = find(~isnan(nan_sub_mat));  %normalized to the number of nontnan components (average)
        
        n=length(notnan_inds);
        
        if n < 20
            out_mat(i,j) = NaN;
            continue;
        end
        
        sigma_x_y =sum(nan_sub_mat(notnan_inds));
        sigma_x =      sum(sub_mat(notnan_inds));
        sigma_y =      sum(mat2(notnan_inds));
        sigma_x2 = sum(sub_mat(notnan_inds).^2);
        sigma_y2 = sum(mat2(notnan_inds).^2);
        
        out_mat(i,j) = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
                                    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
                                    sqrt(n*sigma_y2-sigma_y.^2);
        

         end % for j
    end % for i
disp('')

function out_mat = nanconv2(mat,h)

out_mat = mat;
 nan_mat = isnan(mat);
 
 % dilate nan_mat
 
SE = strel('disk', 2);
nan_mat =  ~imdilate(~nan_mat,SE);

i_size = size(h,1); j_size = size(h,2);
[work_mat,npad_i,npad_j] = pad_edges(mat,h,2);
for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        
        % for each i and j, choose the correct sub-mat (size of h) to multiply with h
        
        sub_mat = work_mat(npad_i+i-floor(i_size/2):npad_i+i+floor(i_size/2), ...
                                                          npad_j+j-floor(j_size/2):npad_j+j+floor(j_size/2)  ); % assumes h is odd in number
                                                      
        notnan_inds = find(~isnan(sub_mat));  
        
        if ~isempty(notnan_inds)
            sum_h = sum(h(notnan_inds));   % normalize to the places without a NaN
            out_mat(i,j) = nansum(nansum(sub_mat .* h));
            out_mat(i,j) = out_mat(i,j)/sum_h;
        end
        
    end % for j
end % for i

out_mat(nan_mat) = NaN;

disp('')
 
function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)

npad_ij = ceil(size(h)/l);
npad_i = npad_ij(1);
npad_j = npad_ij(2);
in_size = size(mat);
out_size = in_size + [2*npad_i 2*npad_j];
out_mat = nan(out_size);
out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;
disp('')
 function score=criteria(sigma)

win_len=ceil(sigma*15);
if sigma==0
win_len=15;
end
[win] = MY_gaussianwin(win_len,sigma);
win=win/sum(win);
%fwvt=wvtool(win);

[val,ind]=max(win);
score=win(ind+1)^2;
if sigma==0
    score=1;
end
disp('')
function [w] = MY_gaussianwin(M,sigma)
n= -(M-1)/2 : (M-1)/2;
w = exp(-n .* n / (2 * sigma * sigma))';
w(isnan(w))=1;
disp('')

function out_mat=PointCorr(Rxx,RxxR)
nan_mat=Rxx .* RxxR;                                             
        notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)
        
        n=length(notnan_inds);
        
        if n < 2
            out_mat = NaN;
            
        end
        
        sigma_x_y =sum(nan_mat(notnan_inds));
        sigma_x =      sum(Rxx(notnan_inds));
        sigma_y =      sum(RxxR(notnan_inds));
        sigma_x2 = sum(Rxx(notnan_inds).^2);
        sigma_y2 = sum(RxxR(notnan_inds).^2);
        
        out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
                                    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
                                    sqrt(n*sigma_y2-sigma_y.^2);
disp('')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xymax,smax,xymin,smin] = extrema2(xy,varargin)
%EXTREMA2   Gets the extrema points from a surface.

%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA2(X) returns the maxima and minima 
%   elements of the matriz X ignoring NaN's, where
%    XMAX - maxima points in descending order (the bigger first and so on)
%    IMAX - linear indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - linear indexes of the XMIN.
%   The program uses EXTREMA.
% 
%   The extrema points are searched only through the column, the row and
%   the diagonals crossing each matrix element, so it is not a perfect
%   mathematical program and for this reason it has an optional argument.
%   The user should be aware of these limitations.
%
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA2(X,1) does the same but without
%   searching through the diagonals (less strict and perhaps the user gets
%   more output points).
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema). 
%
%   Note: To change the linear index to (i,j) use IND2SUB. 
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,3:-.2:-2);
%      z = x.*exp(-x.^2-y.^2); z(10,7)= NaN; z(16:19,13:17) = NaN;
%      surf(x,y,z), shading interp
%      [zmax,imax,zmin,imin] = extrema2(z);
%      hold on  
%       plot3(x(imax),y(imax),zmax,'bo',x(imin),y(imin),zmin,'ro')
%       for i = 1:length(zmax)
%        text(x(imax(i)),y(imax(i)),zmax(i),['  ' num2str(zmax(i))])
%       end
%       for i = 1:length(zmin)
%        text(x(imin(i)),y(imin(i)),zmin(i),['  ' num2str(zmin(i))])
%       end
%      hold off 
%
%   See also EXTREMA, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adrin Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2005
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2006-11-22 : Fixed bug in INDX (by JaeKyu Suhr)
% 2007-04-09 : Change name to MAXIMA2, and definition added.

M = size(xy);
if length(M) ~= 2
 error('Entry must be a matrix.')
end
N = M(2);
M = M(1);

% Search peaks through columns:
[smaxcol,smincol] = extremos(xy);

% Search peaks through rows, on columns with extrema points:
im = unique([smaxcol(:,1);smincol(:,1)]); % Rows with column extrema
[smaxfil,sminfil] = extremos(xy(im,:).');

% Convertion from 2 to 1 index:
smaxcol = sub2ind([M,N],smaxcol(:,1),smaxcol(:,2));
smincol = sub2ind([M,N],smincol(:,1),smincol(:,2));
smaxfil = sub2ind([M,N],im(smaxfil(:,2)),smaxfil(:,1));
sminfil = sub2ind([M,N],im(sminfil(:,2)),sminfil(:,1));

% Peaks in rows and in columns:
smax = intersect(smaxcol,smaxfil);
smin = intersect(smincol,sminfil);

% Search peaks through diagonals?
if nargin==1
 % Check peaks on down-up diagonal:
 [iext,jext] = ind2sub([M,N],unique([smax;smin]));
 [sextmax,sextmin] = extremos_diag(iext,jext,xy,1);

 % Check peaks on up-down diagonal:
 smax = intersect(smax,[M; (N*M-M); sextmax]);
 smin = intersect(smin,[M; (N*M-M); sextmin]);

 % Peaks on up-down diagonals:
 [iext,jext] = ind2sub([M,N],unique([smax;smin]));
 [sextmax,sextmin] = extremos_diag(iext,jext,xy,-1);

 % Peaks on columns, rows and diagonals:
 smax = intersect(smax,[1; N*M; sextmax]);
 smin = intersect(smin,[1; N*M; sextmin]);
end

% Extrema points:
xymax = xy(smax);
xymin = xy(smin);

% Descending order:
[~,inmax] = sort(-xymax); clear temp
xymax = xymax(inmax);
smax = smax(inmax);
[xymin,inmin] = sort(xymin);
smin = smin(inmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [smax,smin] = extremos(matriz)
% Peaks through columns or rows.

smax = [];
smin = [];

for n = 1:length(matriz(1,:))
 [temp,imaxfil,temp,iminfil] = extrema(matriz(:,n)); clear temp
 if ~isempty(imaxfil)     % Maxima indexes
  imaxcol = repmat(n,length(imaxfil),1);
  smax = [smax; imaxfil imaxcol];
 end
 if ~isempty(iminfil)     % Minima indexes
  imincol = repmat(n,length(iminfil),1);
  smin = [smin; iminfil imincol];
 end
end


function [sextmax,sextmin] = extremos_diag(iext,jext,xy,A)
% Peaks through diagonals (down-up A=-1)

[M,N] = size(xy);
if A==-1
 iext = M-iext+1;
end
[iini,jini] = cruce(iext,jext,1,1);
[iini,jini] = ind2sub([M,N],unique(sub2ind([M,N],iini,jini)));
[ifin,jfin] = cruce(iini,jini,M,N);
sextmax = [];
sextmin = [];
for n = 1:length(iini)
 ises = iini(n):ifin(n);
 jses = jini(n):jfin(n);
 if A==-1
  ises = M-ises+1;
 end
 s = sub2ind([M,N],ises,jses);
 [temp,imax,temp,imin] = extrema(xy(s)); clear temp
 sextmax = [sextmax; s(imax)'];
 sextmin = [sextmin; s(imin)'];
end


function [i,j] = cruce(i0,j0,I,J)
% Indexes where the diagonal of the element io,jo crosses the left/superior
% (I=1,J=1) or right/inferior (I=M,J=N) side of an MxN matrix. 

arriba = 2*(I*J==1)-1;

si = (arriba*(j0-J) > arriba*(i0-I));
i = (I - (J+i0-j0)).*si + J+i0-j0;
j = (I+j0-i0-(J)).*si + J;
function [xmax,imax,xmin,imin] = extrema(x)
%EXTREMA   Gets the global extrema points from a time series.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima and maxima 
%   points of the vector X ignoring NaN's, where
%    XMAX - maxima points in descending order
%    IMAX - indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - indexes of the XMIN
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Example:
%      x = 2*pi*linspace(-1,1);
%      y = cos(x) - 0.5 + 0.5*rand(size(x)); y(40:45) = 1.85; y(50:53)=NaN;
%      [ymax,imax,ymin,imin] = extrema(y);
%      plot(x,y,x(imax),ymax,'g.',x(imin),ymin,'r.')
%
%   See also EXTREMA2, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adrin Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2004
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2007-04-09 : Change name to MAXIMA, and definition added.


xmax = [];
imax = [];
xmin = [];
imin = [];

% Vector input?
Nt = numel(x);
if Nt ~= length(x)
 error('Entry must be a vector.')
end

% NaN's:
inan = find(isnan(x));
indx = 1:Nt;
if ~isempty(inan)
 indx(inan) = [];
 x(inan) = [];
 Nt = length(x);
end

% Difference between subsequent elements:
dx = diff(x);

% Is an horizontal line?
if ~any(dx)
 return
end

% Flat peaks? Put the middle element:
a = find(dx~=0);              % Indexes where x changes
lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
d = a(lm) - a(lm-1);          % Number of elements in the flat peak
a(lm) = a(lm) - floor(d/2);   % Save middle elements
a(end+1) = Nt;

% Peaks?
xa  = x(a);             % Serie without flat peaks
b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)  
                        % 0  =>  negative slopes (maxima begin)
xb  = diff(b);          % -1 =>  maxima indexes (but one) 
                        % +1 =>  minima indexes (but one)
imax = find(xb == -1) + 1; % maxima indexes
imin = find(xb == +1) + 1; % minima indexes
imax = a(imax);
imin = a(imin);

nmaxi = length(imax);
nmini = length(imin);                

% Maximum or minumim on a flat peak at the ends?
if (nmaxi==0) && (nmini==0)
 if x(1) > x(Nt)
  xmax = x(1);
  imax = indx(1);
  xmin = x(Nt);
  imin = indx(Nt);
 elseif x(1) < x(Nt)
  xmax = x(Nt);
  imax = indx(Nt);
  xmin = x(1);
  imin = indx(1);
 end
 return
end

% Maximum or minumim at the ends?
if (nmaxi==0) 
 imax(1:2) = [1 Nt];
elseif (nmini==0)
 imin(1:2) = [1 Nt];
else
 if imax(1) < imin(1)
  imin(2:nmini+1) = imin;
  imin(1) = 1;
 else
  imax(2:nmaxi+1) = imax;
  imax(1) = 1;
 end
 if imax(end) > imin(end)
  imin(end+1) = Nt;
 else
  imax(end+1) = Nt;
 end
end
xmax = x(imax);
xmin = x(imin);

% NaN's:
if ~isempty(inan)
 imax = indx(imax);
 imin = indx(imin);
end

% Same size as x:
imax = reshape(imax,size(xmax));
imin = reshape(imin,size(xmin));

% Descending order:
[temp,inmax] = sort(-xmax); clear temp
xmax = xmax(inmax);
imax = imax(inmax);
[xmin,inmin] = sort(xmin);
imin = imin(inmin);


function [x,y,h]=ellipse(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same. 
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same 
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% as a sample of how ellipse works, the following produces a red ellipse
% tipped up at a 45 deg axis from the x axis
% ellipse(1,2,pi/8,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%

% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original 
% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch


% Check the number of input arguments 

if nargin<1,
  ra=[];
end;
if nargin<2,
  rb=[];
end;
if nargin<3,
  ang=[];
end;

%if nargin==1,
%  error('Not enough arguments');
%end;

if nargin<5,
  x0=[];
  y0=[];
end;
 
if nargin<6,
  C=[];
end

if nargin<7,
  Nb=[];
end

% set up the default values

if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;

% work on the variable sizes

x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);

if isstr(C),C=C(:);end;

if length(ra)~=length(rb),
  error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;

% how many inscribed elllipses are plotted

if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end;

% drawing loop

for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    radm=ra(k);
    radn=rb(k);
    if length(ang)==1
      an=ang;
    else
      an=ang(k);
    end;
  elseif length(ra)==1
    xpos=x0(k);
    ypos=y0(k);
    radm=ra;
    radn=rb;
    an=ang;
  elseif length(x0)==length(ra)
    xpos=x0(k);
    ypos=y0(k);
    radm=ra(k);
    radn=rb(k);
    an=ang(k)
  else
    rada=ra(fix((k-1)/size(x0,1))+1);
    radb=rb(fix((k-1)/size(x0,1))+1);
    an=ang(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;

  co=cos(an);
  si=sin(an);
  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
  x=radm*cos(the)*co-si*radn*sin(the)+xpos;
  y=radm*cos(the)*si+co*radn*sin(the)+ypos;
  h(k)=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
  set(h(k),'color',C(rem(k-1,size(C,1))+1,:));

end;


function [semimajor_axis, semiminor_axis, x0, y0, phi] = ellipse_fit(x, y)
%
% ellipse_fit - Given a set of points (x,y), ellipse_fit returns the
% best-fit ellipse (in the Least Squares sense) 
%
% Input:                  
%                       x - a vector of x measurements
%                       y - a vector of y measurements
%
% Output:
%
%                   semimajor_axis - Magnitude of ellipse longer axis
%                   semiminor_axis - Magnitude of ellipse shorter axis
%                   x0 - x coordinate of ellipse center 
%                   y0-  y coordinate of ellipse center 
%                   phi - Angle of rotation in radians with respect to
%                   the x-axis
%
% Algorithm used:
%
% Given the quadratic form of an ellipse: 
%  
%       a*x^2 + 2*b*x*y + c*y^2  + 2*d*x + 2*f*y + g = 0   (1)
%                          
%  we need to find the best (in the Least Square sense) parameters a,b,c,d,f,g. 
%  To transform this into the usual way in which such estimation problems are presented,
%  divide both sides of equation (1) by a and then move x^2 to the
% other side. This gives us:
%
%       2*b'*x*y + c'*y^2  + 2*d'*x + 2*f'*y + g' = -x^2            (2)
%  
%   where the primed parametes are the original ones divided by a.
%  Now the usual estimation technique is used where the problem is
%  presented as:
%
%    M * p = b,  where M = [2*x*y y^2 2*x 2*y ones(size(x))], 
%    p = [b c d e f g], and b = -x^2. We seek the vector p, given by:
%    
%    p = pseudoinverse(M) * b.
%  
%    From here on I used formulas (19) - (24) in Wolfram Mathworld:
%    http://mathworld.wolfram.com/Ellipse.html
%
%
% Programmed by: Tal Hendel <thendel@tx.technion.ac.il>
% Faculty of Biomedical Engineering, Technion- Israel Institute of Technology     
% 12-Dec-2008
%
%--------------------------------------------------------------------------


x = x(:);
y = y(:);

%Construct M
M = [2*x.*y y.^2 2*x 2*y ones(size(x))];

% Multiply (-X.^2) by pseudoinverse(M)
e = M\(-x.^2);

%Extract parameters from vector e
a = 1;
b = e(1);
c = e(2);
d = e(3);
f = e(4);
g = e(5);

%Use Formulas from Mathworld to find semimajor_axis, semiminor_axis, x0, y0
%, and phi

delta = b^2-a*c;

x0 = (c*d - b*f)/delta;
y0 = (a*f - b*d)/delta;

phi = 0.5 * acot((c-a)/(2*b));

% phi = 0.5 * atan((2*b),(c-a));

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));

semimajor_axis = max(a_prime, b_prime);
semiminor_axis = min(a_prime, b_prime);

 if (a_prime < b_prime)
     phi = pi/2 - phi;
 end

disp('')



function R_outer=Find_R_outer(acorr,parms)
    
  % calculate all the extrema points in the spatial autocorrelation 
[zmax,imax,zmin,imin]= extrema2(acorr);
[i,j]=ind2sub(size(acorr),imax);


%put all extrema points in dist
dist(:,1)=j;
dist(:,2)=i;
n=length(i);

%calculate the distance of all extrema to the central peak and put them in
%column 3
dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
 
 % sort the hexonal peaks by distance to the centeral peak
 [score,ind]=sort(dist(:,3));
 dist=dist(ind,:);
 %zmax=zmax(ind);
 R=dist(2,3);
  count=1;
 i=2;
 hex_peaks(1,:,:)=dist(1,:,:);
 
 % finds the first 6 closest peaks to the central peak
 while count<7 && i<=length(dist)
     
     % calculate the min distance of point i from all other already chosen
     % points
     min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
         (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));
     
     % point i needs to be on the right side (we choos only half the point cause its semetrical)
     % and the distance of point i from all other already chosen points
     % needs to be higher than R/2
      if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/1.5)
         hex_peaks(count,:,:)=dist(i,:,:); 
          count=count+1; 
         hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
         hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);
            
          count=count+1; 
        
      end    
       i=i+1;
         
 end
 
  R_outer=max(hex_peaks(:,3))*1.15;  
    
disp('')

function [major_ax,minor_ax,angle,hex_peaks,x0,y0]=Find_Module(acorr,parms)

% calculate all the extrema points in the spatial autocorrelation 
[zmax,imax,zmin,imin]= extrema2(acorr);
[i,j]=ind2sub(size(acorr),imax);


%put all extrema points in dist
dist(:,1)=j;
dist(:,2)=i;
n=length(i);

%calculate the distance of all extrema to the central peak and put them in
%column 3
dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
 
 % sort the hexonal peaks by distance to the centeral peak
 [score,ind]=sort(dist(:,3));
 dist=dist(ind,:);
 %zmax=zmax(ind);
 R=dist(2,3);
  count=1;
 i=2;
 hex_peaks(1,:,:)=dist(1,:,:);
 
 % finds the first 6 closest peaks to the central peak
 while count<7 && i<=length(dist)
     
     % calculate the min distance of point i from all other already chosen
     % points
     min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
         (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));
     
     % point i needs to be on the right side (we choos only half the point cause its semetrical)
     % and the distance of point i from all other already chosen points
     % needs to be higher than R/2
      if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/1.5)
         hex_peaks(count,:,:)=dist(i,:,:); 
          count=count+1; 
         hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
         hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);
            
          count=count+1; 
        
      end    
       i=i+1;
         
 end
   
  % fitting the elipse to the 6 peaks that we have found        
 [major_ax, minor_ax, x0, y0, angle] = ellipse_fit(hex_peaks(1:6,1), hex_peaks(1:6,2));
 
%  %%%%%%%%%%%%%%%%%%%%%%%%%% fixing a bug in ellipse fit %%%%%%%%%%%%%%%%
%     [x1,y1,fig]=ellipse(major_ax,minor_ax,angle,x0,y0,'k',100);
%     [x2,y2,fig]=ellipse(major_ax,minor_ax,-angle,x0,y0,'k',100);
%   
% %     %%%%%%%%%%%%%%%%%%%% calculate which ellipse fits better (distance to points)
% i=1:length(x1);
% d1=0;d2=0;
% for j=1:6
% d1=d1+min((hex_peaks(j,1)-x1(i)).^2+(hex_peaks(j,2)-y1(i)).^2);
% d2=d2+min((hex_peaks(j,1)-x2(i)).^2+(hex_peaks(j,2)-y2(i)).^2);
% end
% %%%%%%%%%%%%%%%%% if d2 smaller put -angel and the bug is fixed
% if d2<d1
%   angle=-angle;  
% end
%  
 
 
 if ~isreal(major_ax) ||  ~isreal(minor_ax)
     major_ax=nan;
     minor_ax=nan;
 end
 
 hex_peaks(7,1)=x0;
 hex_peaks(7,2)=y0;
disp('')