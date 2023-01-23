%GEODESIC solves geodesic problems concerning three types of surface.
%   by using variational method
%   (1) a right cylinder with a circular cross section
%   (2) a right cone with a circular base
%   (3) a sphere
%   The parametric expression and the length of the geodesic can be 
%   obtained. The results are displayed in 3-D plots. Generally, two 
%   curve can be obtained: one is the correct geodesic; the other is
%   the conjugate geodesic which is longer than the correct one.
%
%   r, h:               geometry of the surface.
%   surface:            surface type.
%   theta1, z1, xi1:    the first ending points in local coordinates
%   theta2, z2, xi2:    the second ending points in local coordinates
%   x, y, z:            parametric expression of the geodesic
%   cx,cy,cz:           parametric expression of the conjugate geodesic
%
% Note:
%   1. Symbolic Toolbox is necessary to this program for solving the
%   constants in the parametric expression of conical or spherical
%   geodesic
% Designed by: Lei Wang, <WangLeiBox@hotmail.com>, 29-Nov-2004.
% Last Revision: 10-Dec-2004.
% Dept. Mechanical & Aerospace Engineering, NC State University.
% Copyright (c)2004, Lei Wang <WangLeiBox@hotmail.com>
% $Revision: 1.0 $  $ Date: 10-Dec-2004 10:32:02 $
clear;clc;  close all; warning off
r=1; % Radius
h=2; % Height of a cone
surface = 2; % Surface type
if isempty(ver('symbolic')),
    error('Symbolic Toolbox is MUST. Please make sure it has already been installed.')    
end
figure(1)
switch surface
    case 1,     disp('Geodesic of a Cylindrical surface')
        %%  Geodesic of a Cylindrical surface
        %       x = r*cos(theta); y = r*sin(theta); z = u;
        %       starting points: (theta1,z1); (theta2,z2)
        %------------------------------------------------------------------
        [X,Y,Z] = cylinder(r,500);    
        surf(X,Y,Z);axis equal
        set(findobj('Type','surface'),'FaceColor',[0.85,0.85,0.85],'MeshStyle','row')
        zlim([-0.2 1.1]);
        
        % Set two starting and ending points
        theta1= 0*pi/180;     z1= 0.2; 
        theta2= 90*pi/180;    z2= 0.8; 
        x1= r*cos(theta1);  y1= r*sin(theta1);
        x2= r*cos(theta2);  y2= r*sin(theta2);
        u = linspace(theta1,theta2,500);
 
        % Parametric expression of Geodesic curve for case 1
        x = r*cos(u); y = r*sin(u);
        z = (z2-z1)/(theta2-theta1)*u + (z1*theta2-z2*theta1)/(theta2-theta1);
        Length = sqrt(r^2*(theta2-theta1)^2+(z2-z1)^2);
        % Parametric expression of CONJUGATE Geodesic curve for case 1
        xi1 = theta1+2*pi;  xi2 = theta2; 
        cu = linspace(xi1,xi2,500);
        cx = r*cos(cu); cy = r*sin(cu);        
        cz = (z2-z1)/(xi2-xi1)*cu + (z1*xi2-z2*xi1)/(xi2-xi1);
        
    case 2,     disp('Geodesic of a Conical surface')
        %%  Geodesic of a Conical surface
        %       x = (h-u)/h*r*cos(theta); y = (h-u)/h*r*sin(theta); z = u;
        %       starting points: (theta1,z1); (theta2,z2)
        %------------------------------------------------------------------
        
        [X,Y,Z] = cylinder(-h/r*[0:r/20:r]+h,30); 
        h_cone=surf(X/max(X(:))*r,Y/max(Y(:))*r,Z*h);
        set(h_cone,'FaceColor',[0.85,0.85,0.85],'MeshStyle','both');
        axis equal
        axis([-1.1,1.1,-1.1,1.1,0,2])
        % Set two starting and ending points
        theta1= 0*pi/180;     z1= 0.2;  
        theta2= 90*pi/180;    z2= 1.2; 
        x1= r*cos(theta1)*(h-z1)/h;  y1= r*sin(theta1)*(h-z1)/h;
        x2= r*cos(theta2)*(h-z2)/h;  y2= r*sin(theta2)*(h-z2)/h;
        u = linspace(theta1,theta2,500);
        % Parametric expression of Geodesic curve for case 2
        syms c1 c2
        [c1,c2] = solve(r*(h-z1)*cos(c2-r*theta1/sqrt(r^2+h^2))-c1,...
                        r*(h-z2)*cos(c2-r*theta2/sqrt(r^2+h^2))-c1);
        C1=double(c1(1));  C2=double(c2(1));
        z= h - C1*sec(C2-r*u/sqrt(r^2+h^2))/r;
        x= r*cos(u).*(h-z)/h;  y= r*sin(u).*(h-z)/h;
        Length = C1*h*sqrt(r^2+h^2)/r*(tan(C2-r*theta1/sqrt(r^2+h^2))...
                                      -tan(C2-r*theta2/sqrt(r^2+h^2)))/h^2;
                                  
        % Parametric expression of CONJUGATE Geodesic curve for case 2
        xi1 = theta1+2*pi;  xi2 = theta2; 
        cu = linspace(xi1,xi2,500);
        syms c1 c2        
        [c1,c2] = solve(r*(h-z1)*cos(c2-r*xi1/sqrt(r^2+h^2))-c1,...
                        r*(h-z2)*cos(c2-r*xi2/sqrt(r^2+h^2))-c1);
        C1=double(c1(1));  C2=double(c2(1));
        cz= h - C1*sec(C2-r*cu/sqrt(r^2+h^2))/r;
        cx= r*cos(cu).*(h-cz)/h;  cy= r*sin(cu).*(h-cz)/h;
        
        
    case 3,     disp('Geodesic of a Spherical surface');
        %%  Geodesic of a Spherical surface
        % x = r*sin(phi)*cos(theta); y = r*sin(phi)*sin(theta); z = r*cos(phi);
        %       starting points: (phi1,theta1);(phi1,theta2);
        %------------------------------------------------------------------
        sphere;
        h_surfs=findobj('Type','surface');
        set(h_surfs(1),'FaceColor',[0.85,0.85,0.85],'MeshStyle','both');
        xlim([-1.2 1.2]);ylim([-1.2 1.2]);
        axis equal;
        
        % Set two starting and ending points
        theta1= 0*pi/180;     phi1 = 45*pi/180;
        theta2= 90*pi/180;    phi2 = 100*pi/180;
        x1= r*sin(phi1)*cos(theta1);  y1= r*sin(phi1)*sin(theta1);  z1= r*cos(phi1); 
        x2= r*sin(phi2)*cos(theta2);  y2= r*sin(phi2)*sin(theta2);  z2= r*cos(phi2); 
        u = linspace(theta1,theta2,500); 
        % Parametric expression of Geodesic curve for case 3
        syms c1 c2
        [c1,c2] = solve(sin(c2-theta1)-c1/tan(phi1), sin(c2-theta2)-c1/tan(phi2));
        C1=double(c1(1));  C2=double(c2(1));
        v= atan2(C1,sin(C2-u));
        x= r*sin(v).*cos(u);  y= r*sin(v).*sin(u); z = r*cos(v); i=sqrt(-1);
        Length = r*(atan(sqrt(1+1/C1^2)*tan(C2-theta1))-...
                    atan(sqrt(1+1/C1^2)*tan(C2-theta2)));
        
        % Parametric expression of CONJUGATE Geodesic curve for case 3
        xi1 = theta1+2*pi;  xi2 = theta2; 
        cu = linspace(xi1,xi2,5000);
        syms c1 c2
        [c1,c2] = solve(sin(c2-xi1)-c1/tan(phi1), sin(c2-xi2)-c1/tan(phi2));
        C1=double(c1(1));  C2=double(c2(1));
        cv= atan2(C1,sin(C2-cu));
        cx= r*sin(cv).*cos(cu);  cy= r*sin(cv).*sin(cu); cz = r*cos(cv);
        
        
    otherwise
        error('Please check the input parameter "surface" (1~3)');
end
grid off; hold on;
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('z','fontsize',14);
plot3(x,y,z,'b-','linewidth',2); % Display the geodesic
plot3(cx,cy,cz,'--','linewidth',2,'color',[0,0.5,0]);% Show conjugate geodesic
plot3(x1,y1,z1,'r.','markersize',14); % Show the ending point
plot3(x2,y2,z2,'m.','markersize',14); % Show the ending point
legend('\fontsize{12}Geodesic curve','Conjugate geodesic');
view([150,45]);
%% Numerical length of the geodesic
n_length=0;
for i=1:length(x)-1
    n_length = n_length + norm([x(i+1)-x(i),y(i+1)-y(i),z(i+1)-z(i)],2);
end
disp(['Ending point A(x1,y1,z1): (',num2str(x1),', ',num2str(y1),', ',num2str(z1),')'])
disp(['Ending point B(x2,y2,z2): (',num2str(x2),', ',num2str(y2),', ',num2str(z2),')'])
disp(['   ']);
disp(['Analytical length of the geodesic = ',num2str(Length)])
disp(['Numerical length of the geodesic  = ',num2str(n_length)])
disp(['   ']); disp(datestr(now));