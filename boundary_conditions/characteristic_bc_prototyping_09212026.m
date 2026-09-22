%% prototype characteristic BC (09/21/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

% syms k [3 1] real
% syms v [3 1] real
% syms r p a gm1 real
% 
% rho = r;
% vel = v;
% khat = k;
% half = 1/2;
% 
% V2 = sum(vel.^2);
% xa =1/a;
% xgm1 = 1/gm1;
% 
% b = half*V2*khat + rho*cross(vel,khat);
% H = half*V2 + xgm1*a^2;
% 
% Pmat(1,1) = khat(1);
% Pmat(1,2) = khat(2);
% Pmat(1,3) = khat(3);
% Pmat(1,4) = half*rho*xa;
% Pmat(1,5) = half*rho*xa;
% 
% Pmat(2,1) = vel(1)*khat(1);
% Pmat(2,2) = vel(1)*khat(2) - rho*khat(3);
% Pmat(2,3) = vel(1)*khat(3) + rho*khat(2);
% Pmat(2,4) = half*rho*xa*(vel(1) + khat(1)*a);
% Pmat(2,5) = half*rho*xa*(vel(1) - khat(1)*a);
% 
% Pmat(3,1) = vel(2)*khat(1) + rho*khat(3);
% Pmat(3,2) = vel(2)*khat(2);
% Pmat(3,3) = vel(2)*khat(3) - rho*khat(1);
% Pmat(3,4) = half*rho*xa*(vel(2) + khat(2)*a);
% Pmat(3,5) = half*rho*xa*(vel(2) - khat(2)*a);
% 
% Pmat(4,1) = vel(3)*khat(1) - rho*khat(3);
% Pmat(4,2) = vel(3)*khat(2) + rho*khat(1);
% Pmat(4,3) = vel(3)*khat(3);
% Pmat(4,4) = half*rho*xa*(vel(3) + khat(3)*a);
% Pmat(4,5) = half*rho*xa*(vel(3) - khat(3)*a);
% 
% Pmat(5,1) = b(1);
% Pmat(5,2) = b(2);
% Pmat(5,3) = b(3);
% Pmat(5,4) = half*rho*xa*(H + a*dot(vel,khat));
% Pmat(5,5) = half*rho*xa*(H - a*dot(vel,khat));


syms rho u v w p vn a q g nx ny nz hf real

P(1,1) =  hf*(vn*a + hf*g*q^2)*a;
P(1,2) = -hf*(a*nx + u*g);
P(1,3) = -hf*(a*ny + v*g);
P(1,4) = -hf*(a*nz + w*g);
P(1,5) = hf*g;

P(2,1) =  a^2*nx + a*(w*ny - v*nz) - hf*nx*g*q^2;
P(2,2) =  g*u*nx;
P(2,3) =  a*nz + nx*g*v;
P(2,4) = -a*ny + nx*g*w;
P(2,5) = -g*nx;

P(3,1) =  a^2*ny + a*(u*nz - w*nx) - hf*ny*g*q^2;
P(3,2) = -a*nz + ny*g*u;
P(3,3) =  g*v*ny;
P(3,4) =  a*nx + ny*g*w;
P(3,5) = -g*ny;

P(4,1) =  a^2*nz + a*(v*nx - u*ny) - hf*nz*g*q^2;
P(4,2) =  a*ny + nz*g*u;
P(4,3) = -a*nx + nz*g*v;
P(4,4) =  g*w*nz;
P(4,5) = -g*nz;

P(5,1) =  hf*(-a*vn + hf*g*q^2);
P(5,2) =  hf*( a*nx - g*u );
P(5,3) =  hf*( a*ny - g*v );
P(5,4) =  hf*( a*nz - g*w );
P(5,5) =  hf*g;
