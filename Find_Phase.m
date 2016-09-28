function [brillouin_zone,phase1,Norm_dist,factor]=...
    Find_Phase(hex_peaks_acor,xcor,hex_peaks_xcor)

%finding the Voronoi zone (where to look for the new peak)
brillouin_zone=Find_Voronoi_Zone(xcor,hex_peaks_acor);

%% finding ALL maximum points inside the voronoi zone
[zmax,imax,zmin,imin]= Extrema2(brillouin_zone);
if ~isempty(imax)
    
% finding all local maximum points   
[max_p(:,2),max_p(:,1)]=ind2sub(size(brillouin_zone),imax);

% finding the borders of voronoi zone
% tmp_border = bwboundaries(brillouin_zone);
% Voronoi_border=tmp_border{:,1};
% 
% % subtract all local maximum on the voronoi zone
% %all maximum on the voronoi border
%  max_p_on_borders= intersect(max_p,[Voronoi_border(:,2) Voronoi_border(:,1)],'rows');
%  %% all maximum that are not on the voronoi border
%  [C,i] = setdiff(max_p,max_p_on_borders,'rows');
%  max_p(min(i),1);

C= max_p;

if ~isempty(C)
new_cent(1)=max_p(min(i),1);
 
new_cent(2)=max_p(min(i),2);
 


old_cent(1)=hex_peaks_acor(7,1,1);
old_cent(2)=hex_peaks_acor(7,1,2);


%% line between phase point and center point
x1=new_cent(1);
y1=new_cent(2);

x2=hex_peaks_acor(7,1,1);
y2=hex_peaks_acor(7,1,2);

line_slope=(round(y1)-round(y2))/(round(x1)-round(x2));

if x1-x2>0
    ax=0:0.01:size(xcor,1);
else
    ax=-[0:0.01:size(xcor,1)];
end
%% calculate line
line_y=line_slope*ax+y2;
line_x=ax+x2;

if line_slope==inf
        ax=0:0.01:size(xcor,1);
        line_y=ax+y2;
        line_x=x2*ones(1,length(ax));
end

if line_slope==-inf
        ax=0:0.01:size(xcor,1);
        line_y=-ax+y2;
        line_x=x2*ones(1,length(ax));

end



tmp_border = bwboundaries(brillouin_zone);
Voronoi_border=tmp_border{:,1};

%% calculate the point of intersection between the line and voronoi zone borders
% for i=1:length(line_x)
%     
%  dist_to_line(i)=min(sqrt((line_x(i)-Voronoi_border(:,2)).^2+(line_y(i)-Voronoi_border(:,1)).^2));   
% 
% end
% 
% 
% [val,ind]=min(dist_to_line);

% figure;
% imagesc(brillouin_zone);axis('xy');hold on;
% plot(x2,y2,'*k','linewidth',3);
% plot(new_cent(1),new_cent(2),'*g','linewidth',3);
% plot(line_x,line_y,'*r');
% plot(Voronoi_border(:,2),Voronoi_border(:,1),'*k')
% plot(line_x(ind),line_y(ind),'*b')
% 
% min(line_x)
% max(line_x)
% min(line_y)
% max(line_y)
% 
% 
% norm_factor=sqrt((line_x(ind)-hex_peaks_acor(7,1,1))^2+(line_y(ind)-hex_peaks_acor(7,1,2))^2);
% tmp_phase=sqrt((new_cent(1)-hex_peaks_acor(7,1,1))^2+(new_cent(2)-hex_peaks_acor(7,1,2))^2);
% Norm_dist=tmp_phase/norm_factor*pi;
% 
% 
% factor(1)=line_x(ind);
% factor(2)=line_y(ind);




phase1.score=zmax(1);
phase1.nat_base=new_cent-old_cent;

% trnasform the peak to the module base 
% [phase1.module_base(1),phase1.module_base(2)]=rotatePath...
%     (phase1.nat_base(1),phase1.nat_base(2),-tAngle);
phase1.module_base(3)=sqrt(phase1.module_base(1)^2+phase1.module_base(2)^2);
phase1.dist=sqrt(new_cent(1)^2+new_cent(2)^2);
else % case where all max points are on the border of the Voronoi zone
  phase1.module_base(1)=nan;
  phase1.module_base(2)=nan;
  phase1.nat_base(1)=nan;
  phase1.nat_base(2)=nan;
  phase1.dist=nan;
  Norm_dist=nan;
  factor(1)=nan;
  factor(2)=nan;  
end
else % case where there is no max point inside the Voronoi zone
    phase1.module_base(1)=nan;
    phase1.module_base(2)=nan;
    phase1.nat_base(1)=nan;
    phase1.nat_base(2)=nan;
    phase1.dist=nan;
    Norm_dist=nan;
    factor(1)=nan;
     factor(2)=nan;
end




disp('')



function central_Voronoi_Zone=Find_Voronoi_Zone(xcor,hex_peaks)

%matrix with distance to each peak in the module
size_xcor=size(xcor);
peaks_num= length(hex_peaks);
dist_to_peak=zeros(peaks_num,size_xcor(1),size_xcor(2));


%% calculate the distances matrix from all peaks
[Y,X] = ndgrid(1:1:size(xcor,1), 1:1:size(xcor,2));

for i=1:length(hex_peaks)
    
    
% dist_to_ALL_Peaks(i,:,:)=sqrt((Y-hex_peaks(i,1,2)).^2+(X-hex_peaks(i,1,1)).^2);
 
 dist_to_ALL_Peaks(i,:,:)=sqrt((Y-hex_peaks(i,2)).^2+(X-hex_peaks(i,1)).^2);
 

end


 [dist_to_all_voronoi_z,all_Voronoi_Zones(:,:)]=min(dist_to_ALL_Peaks);
  
 central_Voronoi_Zone=zeros(size_xcor(1),size_xcor(2));
 
 %% voronoi zone 7 is the voronoi of center
 central_Voronoi_Zone(all_Voronoi_Zones==7)=xcor(all_Voronoi_Zones==7);
 
  
disp('') 