%% testing MMS source terms
clc; clear; close all;
load("unsteady_coeffs.mat")

a([8,14,17,20,21,22],:) = 0;
a(:,4) = 0;

% x,y,z,t,in5,l,c,mu,gamma,PrL

PrL = 0.72;
mu  = 0*1.7894e-05;
gamma = 1.4;
c   = 370;
l   = 1.0;

Nx = 100;
Ny = 100;
Nz = 1;

xs = linspace(0,1,Nx);
ys = linspace(0,1,Nx);
zs = linspace(0,1,Nz);

[X,Y,Z] = ndgrid(xs,ys,zs);

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


% levels = [-1.6e07,linspace(-1.4e07,8e06,12)];
% 
% contourf(X,Y,s5(X,Y,Z,0),levels)
% colorbar

contourf(X,Y,uvel(X,Y,Z,0),18)
colorbar