function read_example_haf_layer2(dir_name)
% 
if ~exist('dir_name')
     dir_name = '\\192.114.21.198\Dori_Data\data\rebekkah\Data_sargolini';
end

parms.dir_load_data='\\192.114.21.198\Dori_Data\data\rebekkah\Data_sargolini';
parms.bin_size = 3; % size of spacial bin (for create the rate map)
parms.sigma = 1.5; % gaussiam smoothing factor
parms.time_per_bin=0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
% Minimum radius used in the auto-correlogram when finding the best
% gridness score
parms.minRadius = 20; % [cm]
% Increment step of the radius when calculating what radius gives the best
% gridness score.
parms.radiusStep = 5; % [cm]
parms.pop_grid_size=100; %population size of the shifted grid
parms.interval_num_hist=30;
parms.likelihood_sig=0.1; %threshold for significance gridness likelihood

cell_struct=Read_Examples_2(dir_name,parms);

for i=1:length(cell_struct)
    
    cell_number=i
    
 %   cell_struct(i).rate_mat=Creat_Rate_Map(cell_struct(i).pos.x,cell_struct(i).pos.y,...
 %   cell_struct(i).pos.t,cell_struct(i).spk.x,cell_struct(i).spk.y,cell_struct(i).spk.t,parms);
    
   % cell_struct(i).autocorr=Cross_Correlation(cell_struct(i).rate_mat,cell_struct(i).rate_mat);
 
    % calculate the spacing (ellipse (minor major and phi) and the 6
    % nearest peaks and the center point (x0,y0)
   % [cell_struct(i).spacing.major,cell_struct(i).spacing.minor,cell_struct(i).spacing.phi,...
   %  cell_struct(i).hex_peaks,x0,y0]=Spacing(cell_struct(i).autocorr,parms);
    
    
    
    %finds the autocorrelation for the 6 rotating maps
%     corrMaps = correlationMaps(cell_struct(i).pos.x,cell_struct(i).pos.y,...
%                cell_struct(i).pos.t,cell_struct(i).spk.t,cell_struct(i).rate_mat,parms);
%     
%     % Set the axis for the correlation map
%     corrAxis = parms.bin_size * linspace(-((size(corrMaps{1},1)-1)/2),((size(corrMaps{1},1)-1)/2),size(corrMaps{1},1));

% maxX = nanmax(cell_struct(i).pos.x);
% maxY = nanmax(cell_struct(i).pos.y);
% xStart = nanmin(cell_struct(i).pos.x);
% yStart = nanmin(cell_struct(i).pos.y);
% xLength = maxX - xStart + 10;
% yLength = maxY - yStart + 10;
% tLength = max([xLength,yLength]);
           
% [cell_struct(i).gridness1.value,cell_struct(i).gridness2.value, cell_struct(i).radius1,...
%     cell_struct(i).radius2] = gridnessRadius(corrMaps, tLength, corrAxis,parms);


% [cell_struct(i)]=Gridness_Likelihood(cell_struct(i),parms);
% generate the figure
% [cell_fig]=Ploting_Figure(cell_struct(i),parms);
 

% save the figure and the cells data
S=cell_struct(i);
% if cell_struct(i).gridness1.value>cell_struct(i).gridness1.threshold && cell_struct(i).gridness2.value>cell_struct(i).gridness2.threshold
% cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\grid1grid2\figures\');
% saveas(cell_fig,sprintf('cell_data_%s%s%st%dc%d.jpg',S.rat,S.date,S.session,S.tetrode,S.cell)) ;
 cd('\\192.114.21.198\Dori_Data\data\rebekkah\sargolini data regular format');
 save(sprintf('cell_data_%s%s%st%dc%d.mat',S.rat,S.date,S.session,S.tetrode,S.cell),'S');
 cd(parms.dir_load_data)
end

% if cell_struct(i).gridness1.value>cell_struct(i).gridness1.threshold && cell_struct(i).gridness2.value<cell_struct(i).gridness2.threshold
% % cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\grid1\figures\');
% % saveas(cell_fig,sprintf('cell_data_%s%s%st%dc%d.jpg',S.rat,S.date,S.session,S.tetrode,S.cell)) ;
%  cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\grid1\cells\');
%  save(sprintf('cell_data_%s%s%st%dc%d.mat',S.rat,S.date,S.session,S.tetrode,S.cell),'S');
%  cd(parms.dir_load_data)
% end
% 
% if cell_struct(i).gridness1.value<cell_struct(i).gridness1.threshold && cell_struct(i).gridness2.value>cell_struct(i).gridness2.threshold
% % cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\grid2\figures\');
% % saveas(cell_fig,sprintf('cell_data_%s%s%st%dc%d.jpg',S.rat,S.date,S.session,S.tetrode,S.cell)) ;
%  cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\grid2\cells\');
%  save(sprintf('cell_data_%s%s%st%dc%d.mat',S.rat,S.date,S.session,S.tetrode,S.cell),'S');
%  cd(parms.dir_load_data)
% end
% 
% if cell_struct(i).gridness1.value<cell_struct(i).gridness1.threshold && cell_struct(i).gridness2.value<cell_struct(i).gridness2.threshold
% % cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\no grid1 no grid2\figures\');
% % saveas(cell_fig,sprintf('cell_data_%s%s%st%dc%d.jpg',S.rat,S.date,S.session,S.tetrode,S.cell)) ;
%  cd('C:\Users\Dori\Desktop\raw data\gridcell_data_sargolini\data base\no grid1 no grid2\cells\');
%  save(sprintf('cell_data_%s%s%st%dc%d.mat',S.rat,S.date,S.session,S.tetrode,S.cell),'S');
% cd(parms.dir_load_data)
% end
% end




disp('')


function [cell_struct]=Gridness_Likelihood(cell_struct,parms)

n=parms.pop_grid_size;
for i=1:n
    
    new_spkt=Shift_Spike_Times(cell_struct.spk.t,cell_struct.pos.t);
    
    % adujst the spike posiotion so they will fit the rat path
    spkx=interp1(cell_struct.pos.t,cell_struct.pos.x,new_spkt);
    spky=interp1(cell_struct.pos.t,cell_struct.pos.y,new_spkt);
    
    
    % generate new rate map
    map = Creat_Rate_Map(cell_struct.pos.x,cell_struct.pos.y,cell_struct.pos.t,...
        spkx,spky,new_spkt,parms);
    
    %finds the autocorrelation for the 6 rotating maps
    corrMaps = correlationMaps(cell_struct.pos.x,cell_struct.pos.y,...
        cell_struct.pos.t,new_spkt,map,parms);
    
    % Set the axis for the correlation map
    corrAxis = parms.bin_size * linspace(-((size(corrMaps{1},1)-1)/2),((size(corrMaps{1},1)-1)/2),size(corrMaps{1},1));

    maxX = nanmax(cell_struct.pos.x);
    maxY = nanmax(cell_struct.pos.y);
    xStart = nanmin(cell_struct.pos.x);
    yStart = nanmin(cell_struct.pos.y);
    xLength = maxX - xStart + 10;
    yLength = maxY - yStart + 10;
    tLength = max([xLength,yLength]);
    [gridness1(i),gridness2(i),radius1,radius2] = gridnessRadius(corrMaps, tLength, corrAxis,parms);
     
    i
end
cell_struct.gridness1.pop=gridness1;
min_grid=min(gridness1);
max_grid=max(gridness1);
lag=(max_grid-min_grid)/(parms.interval_num_hist-1);
ax=min_grid:lag:max_grid;
h=hist(gridness1,ax);
h=h/sum(h);
cell_struct.gridness1.hist1=h;
cell_struct.gridness1.hist1ax=ax;

i=1;
while(sum(h(1:i))<1-parms.likelihood_sig)
    i=i+1;
end
cell_struct.gridness1.threshold=ax(i);

cell_struct.gridness2.pop=gridness2;
min_grid=min(gridness2);
max_grid=max(gridness2);
lag=(max_grid-min_grid)/(parms.interval_num_hist-1);
ax=min_grid:lag:max_grid;
h=hist(gridness2,ax);
h=h/sum(h);
cell_struct.gridness2.hist2=h;
cell_struct.gridness2.hist2ax=ax;
i=1;
while(sum(h(1:i))<1-parms.likelihood_sig)
    i=i+1;
end
cell_struct.gridness2.threshold=ax(i);

 
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

% Calculates the correlogram for the different rotations
function corrMaps = correlationMaps(posx,posy,post,spkt,map,parms)

% Number of rotations 180/30=6
rotStep = 30;
N = floor(180/rotStep);

corrMaps = cell(N,1);
corrMaps{1} =Cross_Correlation(map, map);

% Number of bins in the correlogram
Nr = size(corrMaps{1},1);


rotAngle = 0;
for ii = 1:N-1

%     figure;imagesc(map);title('rate map');axis equal
    
      
    % Increment rotation angle
    rotAngle = rotAngle + rotStep;

    % Convert angle from degrees to radians
    radAngle = rotAngle*pi/180;
    
    % rotate the rat path (posx,poy)
    [rx,ry] = rotatePath(posx,posy,radAngle);
    
    % adujst the spike posiotion so they will fit the rat path
    spkx=interp1(post,rx,spkt);
    spky=interp1(post,ry,spkt);
    
   
    %creat the new rate map
    map = Creat_Rate_Map(rx,ry,post,spkx,spky,spkt,parms);
    

    
    % Calculate the autocorrelation map for the rotated data
    corrMaps{ii+1} = Cross_Correlation(map,map);
    
    corrMaps{ii+1} = adjustMapSize(corrMaps{ii+1},Nr);

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rxx_new = adjustMapSize(Rxx,Nr)

% Size of map, new size must be Nr
N = size(Rxx,1);

 if N == Nr 
     Rxx_new = Rxx;
     return
 end

Rxx_new = nan(Nr);

if N > Nr
    diff = (N - Nr) / 2;
    for ii = 1:Nr
        
        for jj = 1:Nr
            Rxx_new(ii,jj) = Rxx(ii+diff,jj+diff);
           
        end
    end
    
end

if N < Nr
    diff = (Nr - N) / 2;
    for ii = 1:N
        for jj = 1:N
            Rxx_new(ii+diff,jj+diff) = Rxx(ii,jj);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated the position angles according to the angle tAngle [radians]
function [newX,newY] = rotatePath(posx,posy,tAngle)

% multiply by a rotation matrix
newX = posx * cos(tAngle) - posy * sin(tAngle); %rotate posx
newY = posx * sin(tAngle) + posy * cos(tAngle); % rotate posy


% Calculates the radius that gives the best gridness score
function [gridness1,gridness2, radius1,radius2] = gridnessRadius(corrMaps, tLength, corrAxis,parms)


N = floor(180/30);

numRadius = floor((tLength-parms.minRadius)/parms.radiusStep);

% Array to hold the gridness scores for the different radiuses
gridnessArray1 = zeros(numRadius,2);
gridnessArray2 = zeros(numRadius,2);
corrValues = zeros(N,numRadius);

% Set the radius
radius = parms.minRadius;

for ii = 1:numRadius
    
    
    Rxx = corrMaps{1};
    % Set the part of the map outside the radius to NaN
    Rxx = adjustMap(Rxx,radius,-1,corrAxis);
    
    for jj = 2:N
        RxxR = corrMaps{jj};
        
        % Set the part of the map outside the radius to NaN
        RxxR = adjustMap(RxxR,radius,-1,corrAxis);
        
        corrValues(jj,ii) = PointCorr(Rxx,RxxR);

    end
    
    % Calculate the degree of "gridness"
    
    % difference of means
    sTop = mean([corrValues(3,ii),corrValues(5,ii)]);
    sTrough = mean([corrValues(2,ii),corrValues(4,ii),corrValues(6,ii)]);
    gridnessArray1(ii,1) = sTop - sTrough;

    % below definition is done according to langston et al. (2010)
    sTop_vec = [corrValues(3,ii),corrValues(5,ii)];
    sTrough_vec = [corrValues(2,ii),corrValues(4,ii),corrValues(6,ii)];
    diff_mat = repmat(sTop_vec',1,3)-repmat(sTrough_vec,2,1);
    gridnessArray2(ii,1) = min(diff_mat(:));
    
    gridnessArray1(ii,2) = radius;
    gridnessArray2(ii,2) = radius;


    
    % Increment the radius
    radius = radius + parms.radiusStep;
end


% Maximum gridness
[gridness1,gInd] = max(gridnessArray1(:,1));
radius1 = gridnessArray1(gInd,2);
[gridness2,gInd] = max(gridnessArray2(:,1));
radius2 = gridnessArray2(gInd,2);
corrValues = corrValues(:,gInd);

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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets the bins of the map outside the radius to NaN
function Rxx = adjustMap(Rxx,radius,centreRadius,mapAxis)

% Size, in number of bins, of the map
[N,M] = size(Rxx);


% Calculate the distance from origo for each bin in the map
oDist = zeros(N,M);
for ii = 1:N
    for jj = 1:M
        oDist(ii,jj) = sqrt(mapAxis(ii)^2 + mapAxis(jj)^2);
       
    end
   
end
Rxx(oDist>radius) = NaN;
Rxx(oDist<=centreRadius) = NaN;



function [major_ax,minor_ax,angle,hex_peaks,x0,y0]=Spacing(acorr,parms)


%  h=fspecial('gaussian',5*[parms.sigma,parms.sigma],0.5*parms.sigma);
%   
%  %smooth the rate mat
% smooth_acorr=nanconv2(acorr,h);
%find the local maxima in the autocorr matrix
[zmax,imax,zmin,imin]= extrema2(acorr);
[i,j]=ind2sub(size(acorr),imax);

%calculat the distance from the central peak to the other peaks
dist(:,1)=j;
dist(:,2)=i;
n=length(i);
dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
 
 % sort and the hexonal peaks by distance to the center
 [score,ind]=sort(dist(:,3));
 dist=dist(ind,:);
 R=dist(2,3);
  count=1;
 i=2;
 hex_peaks(1,:,:)=dist(1,:,:);
 while count<7
     min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
         (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));
     
      if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/2)
         hex_peaks(count,:,:)=dist(i,:,:); 
          count=count+1; 
         hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
         hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);
         i      
          count=count+1; 
        
      end    
       i=i+1;
         
 end
   
          
 [major_ax, minor_ax, x0, y0, angle] = ellipse_fit(hex_peaks(1:6,1), hex_peaks(1:6,2));
 hex_peaks(7,1)=x0;
 hex_peaks(7,2)=y0;
 
function [h]=Ploting_Figure(cell_struct,parms)

 set(0,'defaulttextinterpreter','none');
 


 m=figure;
 s2= sprintf('dir: %s\nfile: %s',cell_struct.dir_name,cell_struct.file_name);
 subplot(2,3,1);plot(cell_struct.pos.x,cell_struct.pos.y,...
   'k',cell_struct.spk.x,cell_struct.spk.y,'.R');hold on;title(s2) ;axis equal;axis off
 subplot(2,3,2);imagesc(cell_struct.rate_mat);axis('xy');axis equal;axis off; hold on;
 id=sprintf('Date:%s\nRat:%s\nSession:%s\nTetrode: %d\nCell: %d',...
    cell_struct.date,cell_struct.rat,cell_struct.session,cell_struct.tetrode,cell_struct.cell);
title(id);
 
subplot(2,3,3);
h=imagesc(cell_struct.autocorr);axis('xy'); axis equal;axis off;hold on;
s1= sprintf('Module:\nA=%f\nB=%f\nphi=%f',cell_struct.spacing.major,cell_struct.spacing.minor,cell_struct.spacing.phi);
title(s1);
x0=cell_struct.hex_peaks(7,1);
y0=cell_struct.hex_peaks(7,2);
 beta = cell_struct.spacing.phi ;
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha =0:pi/100:2*pi;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    x = x0 + (cell_struct.spacing.major * cosalpha * cosbeta - cell_struct.spacing.minor * sinalpha * sinbeta);
    y = y0 + (cell_struct.spacing.major * cosalpha * sinbeta + cell_struct.spacing.minor * sinalpha * cosbeta);

 plot(x,y,'k','LineWidth',1);hold on;
 plot(cell_struct.hex_peaks(:,1), cell_struct.hex_peaks(:,2),'ok','LineWidth',1);hold on;
 
 xMajor1 = x0 + cell_struct.spacing.major * cos(cell_struct.spacing.phi);
xMajor2 = x0 - cell_struct.spacing.major * cos(cell_struct.spacing.phi);
yMajor1 = y0 + cell_struct.spacing.major * sin(cell_struct.spacing.phi);
yMajor2 = y0 - cell_struct.spacing.major * sin(cell_struct.spacing.phi);
 p1=xMajor1:(xMajor2-xMajor1)/10:xMajor2;
  p2=yMajor1:(yMajor2-yMajor1)/10:yMajor2;
plot(p1,p2,'LineWidth',1); hold on;

 xMinor1 = x0 + cell_struct.spacing.minor * cos(cell_struct.spacing.phi+pi/2);
xMinor2 = x0 - cell_struct.spacing.minor * cos(cell_struct.spacing.phi+pi/2);
yMinor1 = y0 + cell_struct.spacing.minor * sin(cell_struct.spacing.phi+pi/2);
yMinor2 = y0 - cell_struct.spacing.minor * sin(cell_struct.spacing.phi+pi/2);

p11=xMinor1:(xMinor2-xMinor1)/10:xMinor2;
p21=yMinor1:(yMinor2-yMinor1)/10:yMinor2;
plot(p11,p21,'k','LineWidth',1); hold on;

 i=1;
 precent=0
ax_val=cell_struct.gridness1.hist1ax(i);
while ax_val<cell_struct.gridness1.value && i<=length(cell_struct.gridness1.hist1ax)
    precent=sum(cell_struct.gridness1.hist1(1:i));
    ax_val=cell_struct.gridness1.hist1ax(i);
    i=i+1;
end
subplot(2,3,4);bar(cell_struct.gridness1.hist1ax,cell_struct.gridness1.hist1,'k');hold on;
g1=sprintf('Gridness1=%f (L=%g) , sig(%g)=%f',cell_struct.gridness1.value,1-precent,parms.likelihood_sig,cell_struct.gridness1.threshold);
title(g1);
m=0:0.001:0.07;
plot(cell_struct.gridness1.threshold,m,'.r');hold on
plot(cell_struct.gridness1.value,m,'.g');

i=1;
precent=0;
ax_val=cell_struct.gridness2.hist2ax(i);
while ax_val<cell_struct.gridness2.value && i<=length(cell_struct.gridness2.hist2ax)
    precent=sum(cell_struct.gridness2.hist2(1:i));
    ax_val=cell_struct.gridness2.hist2ax(i);
    i=i+1;
end
subplot(2,3,6);bar(cell_struct.gridness2.hist2ax,cell_struct.gridness2.hist2,'k');hold on;
g2=sprintf('Gridness2=%f(L=%g) , sig(%g)=%f',cell_struct.gridness2.value,...
    1-precent,parms.likelihood_sig,cell_struct.gridness2.threshold);
title(g2);
plot(cell_struct.gridness1.threshold,m,'.r');hold on
 plot(cell_struct.gridness1.value,m,'.g');
  
 




 

 
 

 

 
 



function cell_struct=Read_Examples_2(dir_name,parms)
dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};



% dir_list = dir(dir_name);
% 
% file_names = {dir_list.name};

nfiles = length(file_names);
 
% read cell structure

cell_ind = 0;

for nfile = 1:nfiles
   
    file_name = file_names{nfile};
    
    % parse file name
    
    % determine file type
 
    if (length(file_name) > 8 && file_name(end-7) == 'T' && file_name(end-5) == 'C')
        % extract cell parameters
        
        rat = file_name(1:5);
        date = file_name(7:12);
        suffix = file_name(13:end);
        ind = find(suffix == '_');
        session = suffix(1:ind-1);
        tetrode = str2num(suffix(ind+2));
        cell = str2num(suffix(ind+4));
                spike_data = load([dir_name '/' file_name]);

        
        cell_ind = cell_ind+1;
        
        cell_struct(cell_ind).rat = rat;
        cell_struct(cell_ind).date = date;
        cell_struct(cell_ind).session = session;
        cell_struct(cell_ind).tetrode = tetrode;
        cell_struct(cell_ind).cell = cell;
        cell_struct(cell_ind).spk.t = spike_data.cellTS;
        
        % load correspoding pos file
        ind_=find(file_name=='_');
        pos_file_name = sprintf('%s_POS.mat',file_name(1:ind_-1));
        pos_data = load(pos_file_name);
        cell_struct(cell_ind).pos.x =pos_data.posx;
        cell_struct(cell_ind).pos.y =pos_data.posy;
        
        if isfield(pos_data,'posx2') && isfield(pos_data,'posy2')
            cell_struct(cell_ind).pos.x2 =pos_data.posx2;
            cell_struct(cell_ind).pos.y2 =pos_data.posy2;
        end
            
            
        cell_struct(cell_ind).pos.t =pos_data.post;
        cell_struct(cell_ind).dir_name='gridcell_data_sargolini';
        cell_struct(cell_ind).file_name=file_name;
        %%%%%%%
        %update parms.time_per_bin
        parms.time_per_bin = cell_struct(cell_ind).pos.t(2)-cell_struct(cell_ind).pos.t(1); % in seconds
        
        % interpolation of spik times to possition of the rat
        cell_struct(cell_ind).spk.x=interp1(cell_struct(cell_ind).pos.t,cell_struct(cell_ind).pos.x,cell_struct(cell_ind).spk.t);
        cell_struct(cell_ind).spk.y=interp1(cell_struct(cell_ind).pos.t,cell_struct(cell_ind).pos.y,cell_struct(cell_ind).spk.t);
        
%        figure;
%         plot(cell_struct(cell_ind).pos.x,cell_struct(cell_ind).pos.y,'k',...
%             cell_struct(cell_ind).spk.x,cell_struct(cell_ind).spk.y,'R.');
%         legend('trajectories','spikes');xlabel('X (cm');ylabel('y (cm)');
       
    end
    
end


disp('')

function cell_struct=Creat_ALL_Rate_Maps(cell_struct,parms)
for i=1:length(cell_struct)
    cell_struct(i).rate_mat=Creat_Rate_Map(cell_struct(i),parms);
    i
end

function rate_mat=Creat_Rate_Map(posx,posy,post,spkx,spky,spkt,parms)
max_x = ceil(max(posx)); 
max_y = ceil(max(posy));
min_x = min(floor(posx));
min_y = min(floor(posy));

% divid the environment into spatial bins 
axis_x = min_x:parms.bin_size:max_x;
axis_y = min_y:parms.bin_size:max_y;

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
        time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind)+parms.time_per_bin;
        
    end
end
%create conut mat( count the num of spikes in each bin)
for i = 1:length(spkt)
   
        [min_val,x_ind] =  min(abs(spkx(i)-axis_x));
        [min_val,y_ind] =  min(abs(spky(i)-axis_y));
        spike_mat(y_ind,x_ind)= spike_mat(y_ind,x_ind)+1;
end

% create rate mat
 rate_mat=spike_mat./time_mat;
 rate_mat(rate_mat==inf)=NaN;
 rate_mat=Smooth_Rate_Mat(rate_mat,parms);

function rate_mat=Smooth_Rate_Mat(rate_mat,parms) 
 %create window
  h=fspecial('gaussian',3*[parms.sigma,parms.sigma],parms.sigma);
  
 %smooth the rate mat
rate_mat=nanconv2(rate_mat,h);
%figure;imagesc(rate_mat);title('rate map');axis equal

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

function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)

npad_ij = ceil(size(h)/l);
npad_i = npad_ij(1);
npad_j = npad_ij(2);
in_size = size(mat);
out_size = in_size + [2*npad_i 2*npad_j];
out_mat = nan(out_size);
out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;

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
[temp,inmax] = sort(-xymax); clear temp
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

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));

semimajor_axis = max(a_prime, b_prime);
semiminor_axis = min(a_prime, b_prime);

if (a_prime < b_prime)
    phi = pi/2 - phi;
end