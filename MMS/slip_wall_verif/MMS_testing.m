%% testing MMS source terms (07/29/2025)
clc; clear; close all;
load("unsteady_coeffs.mat")

% a([8,14,17,20,21,22],:) = 0;
% a(:,4) = 0;
a(:,4) = a(:,2) + a(:,3);

% x,y,z,t,in5,l,c,mu,gamma,PrL

PrL = 0.72;
mu  = 0*1.7894e-05;
gamma = 1.4;
c   = 370;
l   = 1.0;

Nx = 129;
Ny = 129;
Nz = 1;

% xs = linspace(0,1,Nx);
% ys = linspace(0,1,Nx);
% zs = linspace(0,1,Nz);
% [X,Y,Z] = ndgrid(xs,ys,zs);
ang = pi/8;
th1 = linspace(-ang,ang,Nx) + pi/2;
th2 = linspace(-ang,ang,Ny) + pi/2;
r   = linspace(0.5,1.1,Nz);
[TH1,TH2,R] = ndgrid(th1,th2,r);

X = R.*sin(TH1).*cos(TH2);
Y = R.*sin(TH1).*sin(TH2);
Z = R.*cos(TH1);

% rho  = @(x,y,z,t) arrayfun(@(x,y,z) dens_mms(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% uvel = @(x,y,z,t) arrayfun(@(x,y,z) uvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% vvel = @(x,y,z,t) arrayfun(@(x,y,z) vvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% wvel = @(x,y,z,t) arrayfun(@(x,y,z) wvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% p    = @(x,y,z,t) arrayfun(@(x,y,z) pres_mms(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
rho  = @(x,y,z,t) dens_mms(x,y,z,t,a,l,c,mu,gamma,PrL);
uvel = @(x,y,z,t) uvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL);
vvel = @(x,y,z,t) vvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL);
wvel = @(x,y,z,t) wvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL);
p    = @(x,y,z,t) pres_mms(x,y,z,t,a,l,c,mu,gamma,PrL);

vmag = @(x,y,z,t) sqrt(  uvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL).^2 ...
                       + vvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL).^2 ...
                       + wvel_mms(x,y,z,t,a,l,c,mu,gamma,PrL).^2 );

% s1 = @(x,y,z,t) arrayfun(@(x,y,z) mass_source(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% s2 = @(x,y,z,t) arrayfun(@(x,y,z) xmtm_source(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% s3 = @(x,y,z,t) arrayfun(@(x,y,z) ymtm_source(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% s4 = @(x,y,z,t) arrayfun(@(x,y,z) zmtm_source(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );
% s5 = @(x,y,z,t) arrayfun(@(x,y,z) ener_source(x,y,z,t,a,l,c,mu,gamma,PrL), x, y, z );

s1 = @(x,y,z,t) mass_source(x,y,z,t,a,l,c,mu,gamma,PrL);
s2 = @(x,y,z,t) xmtm_source(x,y,z,t,a,l,c,mu,gamma,PrL);
s3 = @(x,y,z,t) ymtm_source(x,y,z,t,a,l,c,mu,gamma,PrL);
s4 = @(x,y,z,t) zmtm_source(x,y,z,t,a,l,c,mu,gamma,PrL);
s5 = @(x,y,z,t) ener_source(x,y,z,t,a,l,c,mu,gamma,PrL);

vn = @(x,y,z,t) v_dot_n(x,y,z,t,a,l,c,mu,gamma,PrL);


% levels = [-1.6e07,linspace(-1.4e07,8e06,12)];
% 
% contourf(X,Y,s5(X,Y,Z,0),levels)
% colorbar

% contourf(X,Y,uvel(X,Y,Z,0),18)
% surf(X,Y,Z,uvel(X,Y,Z,0),'EdgeColor','none')
surf(squeeze(X),squeeze(Y),squeeze(Z),squeeze(vn(X,Y,Z,0)),'EdgeColor','none')
colorbar
hold on
quiver3(X,Y,Z,uvel(X,Y,Z,0),vvel(X,Y,Z,0),wvel(X,Y,Z,0),2,"filled",'k')
% quiver3(X,Y,Z,s2(X,Y,Z,0),s3(X,Y,Z,0),s4(X,Y,Z,0))
axis equal