%% testing MMS source terms 08032026
clc; clear; close all;

PrL = 0.72;
mu  = 0*1.7894e-05;
gamma = 1.4;
c   = 370;
l   = 1.0;

Nx = 100;
Ny = 100;
Nz = 1;

xs_ = linspace(0.5,1.0,Nx);
ys_ = linspace(0.0,0.5,Nx);
zs_ = linspace(0,1,Nz);

[X_,Y_,Z] = ndgrid(xs_,ys_,zs_);
X = (4/3)*(X_.^2-(1/4))+1;
Y = Y_ + 0.05*sin(2*pi*X);

rho  = @(x,y,z,t) dens_mms(x,y,z,t,mu,gamma,PrL);
uvel = @(x,y,z,t) uvel_mms(x,y,z,t,mu,gamma,PrL);
vvel = @(x,y,z,t) vvel_mms(x,y,z,t,mu,gamma,PrL);
wvel = @(x,y,z,t) wvel_mms(x,y,z,t,mu,gamma,PrL);
p    = @(x,y,z,t) pres_mms(x,y,z,t,mu,gamma,PrL);

s1 = @(x,y,z,t) mass_source(x,y,z,t,mu,gamma,PrL);
s2 = @(x,y,z,t) xmtm_source(x,y,z,t,mu,gamma,PrL);
s3 = @(x,y,z,t) ymtm_source(x,y,z,t,mu,gamma,PrL);
s4 = @(x,y,z,t) zmtm_source(x,y,z,t,mu,gamma,PrL);
s5 = @(x,y,z,t) ener_source(x,y,z,t,mu,gamma,PrL);

% contourf(X,Y,s5(X,Y,Z,0))
% colorbar

contourf(X,Y,sqrt((uvel(X,Y,Z,0),18)
colorbar