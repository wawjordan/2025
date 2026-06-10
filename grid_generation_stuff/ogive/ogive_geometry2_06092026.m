%% ogive geometry work (06/09/2026) (no symbolic shit)
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

rho_ogv = 2.8;
L0      = 0.5;
Rn      = 12.5/1000;
Rb      = 45/1000;
y_ogive = @(x) sqrt(rho_ogv^2 - (L0-x).^2) + Rb - rho_ogv;
dydx_ogive = @(x) (L0-x)./sqrt(rho_ogv^2 - (L0 - x).^2);
y_circ    = @(x,x0) sqrt( Rn^2 - (x-x0).^2 );
dydx_circ = @(x,x0) -(x-x0)./sqrt(Rn^2- (x-x0).^2);


hold on
fplot(y_ogive,[0,0.1],'r')
axis equal

x01 = L0 + sqrt( (rho_ogv-Rn)^2 - (Rb-rho_ogv)^2 );
x02 = L0 - sqrt( (rho_ogv-Rn)^2 - (Rb-rho_ogv)^2 );

x0_guess = x02;
fplot(@(x)y_circ(x,x0_guess),[x0_guess-Rn,x0_guess+Rn],'k')

x0_guess = 0.075;
fplot(@(x)y_circ(x,x0_guess),[x0_guess-Rn,x0_guess+Rn],'b')


% 
