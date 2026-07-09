%% orthogonal projection onto airfoil (07/08/2026)
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

airfoil = local_airfoil_generator();

airfoil.plot_airfoil(true);

x_in = 0.8;
y_in = 0.01;

plot(x_in,y_in,'ro');
theta_guess = pi;
[theta,x_out,y_out] = airfoil.point_orthprj_theta(x_in,y_in,2*pi-0.5,true);

plot(x_out,y_out,'bx')


function airfoil = local_airfoil_generator()
epsilon = 0.1;
kappa   = 0.1;
tau     = 0.0;
vinf    = 75.0;
rhoinf  = 1.0;
pinf    = 100000.0;
gamma   = 1.4;
alpha   = 0; % (degrees)
rho_ref = 1.0;
p_ref   = 100000.0;
a_ref   = sqrt(gamma*p_ref/rho_ref);

% nondimensionalize inputs
vinf   = vinf/a_ref;
rhoinf = rhoinf/rho_ref;
pinf   = pinf/(rho_ref*a_ref^2);

airfoil        = kt_airfoil( epsilon, kappa, tau );
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(alpha);
end