function [fig_pair]= plot_ellipse(rate_mat_after_auto_pearson, module_major, module_minor ,module_phi,six_orientation_pts, fig_pair)
 

%mod2= sprintf('A=%0.3f\nB=%0.3f\nphi=\n%0.3f(rad)\n%0.3f(deg)',major,minor,phi,radtodeg(phi));

 lw=2;
 imagesc(rate_mat_after_auto_pearson);axis('xy');axis equal;axis off;hold on;
 axis ij;
 
%center point
x0=six_orientation_pts(7,1); 
y0=six_orientation_pts(7,2);
 beta =module_phi ;
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha =0:pi/100:2*pi;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha); 

    x1 = x0 + (module_major * cosalpha * cosbeta - module_minor * sinalpha * sinbeta);
    y1 = y0 + (module_major * cosalpha * sinbeta + module_minor * sinalpha * cosbeta);

 plot(x1,y1,'k','LineWidth',lw);hold on;
 plot(six_orientation_pts(:,2),six_orientation_pts(:,1),'ok','LineWidth',lw);hold on;
 
xMajor1 = x0 + module_major * cos(module_phi);
xMajor2 = x0 - module_major * cos(module_phi);
yMajor1 = y0 + module_major * sin(module_phi);
yMajor2 = y0 - module_major * sin(module_phi);
 p1=xMajor1:(xMajor2-xMajor1)/10:xMajor2;
 p2=yMajor1:(yMajor2-yMajor1)/10:yMajor2;
plot(p1,p2,'LineWidth',lw); hold on;

xMinor1 = x0 + module_minor * cos(module_phi+pi/2);
xMinor2 = x0 - module_minor * cos(module_phi+pi/2);
yMinor1 = y0 + module_minor * sin(module_phi+pi/2);
yMinor2 = y0 - module_minor * sin(module_phi+pi/2);

p11=xMinor1:(xMinor2-xMinor1)/10:xMinor2;
p21=yMinor1:(yMinor2-yMinor1)/10:yMinor2;
plot(p11,p21,'k','LineWidth',lw); hold on;

disp('')