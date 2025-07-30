%% testing SVF solution (07/30/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma  = 1.4;
inputs = [2.0, 1.0, 0.8611583247416177E+05, 0.5];
% inputs = [2.0, 1.0, 0.8611583247416177E+05, 2.0;
ref_inputs = [1.0, 1.0, 347.2206293753677000 ];

N_t = 33;
N_r = 33;
N_z = 1;

t_range = [pi/2, 0];
r_range = [2, 3];
z_range = [0, 0];

t = linspace(t_range(1),t_range(2),N_t);
r = linspace(r_range(1),r_range(2),N_r);
z = linspace(z_range(1),z_range(2),N_z);
[T,R,Z] = ndgrid(t,r,z);
X = R.*cos(T);
Y = R.*sin(T);

rho  = @(x,y,z) dens_svf(x,y,z,gamma,inputs(:),ref_inputs(:));
uvel = @(x,y,z) uvel_svf(x,y,z,gamma,inputs(:),ref_inputs(:));
vvel = @(x,y,z) vvel_svf(x,y,z,gamma,inputs(:),ref_inputs(:));
wvel = @(x,y,z) wvel_svf(x,y,z,gamma,inputs(:),ref_inputs(:));
p    = @(x,y,z) pres_svf(x,y,z,gamma,inputs(:),ref_inputs(:));

vmag = @(x,y,z) sqrt(  uvel_svf(x,y,z,gamma,inputs(:),ref_inputs(:)).^2 ...
                     + vvel_svf(x,y,z,gamma,inputs(:),ref_inputs(:)).^2 ...
                     + wvel_svf(x,y,z,gamma,inputs(:),ref_inputs(:)).^2 );

contourf(X,Y,rho(X,Y,Z),18)
% surf(X,Y,Z,uvel(X,Y,Z,0),'EdgeColor','none')
% surf(squeeze(X),squeeze(Y),squeeze(Z),squeeze(rho(X,Y,Z)),'EdgeColor','none')
colorbar
hold on
% quiver3(X,Y,Z,uvel(X,Y,Z),vvel(X,Y,Z),wvel(X,Y,Z),2,"filled",'k')
% axis equal