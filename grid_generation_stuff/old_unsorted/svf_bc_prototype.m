%% Prototyping for exact solution face states at boundary (04/26/2023)
clc; clear; close all;

% this --> mms_t (svf)
% grid --> assume center is at (0,0)????
% this%radius_inner

% GRID.x = R.*cos(T);
% GRID.y = R.*sin(T);


syms rho_i gamma M_i r_i r

rho(r) = rho_i * (1 + (7/5-1)/2 * M_i^2 *(1 - r_i^2/r^2))^(1/(7/5-1));

int(rho,r)