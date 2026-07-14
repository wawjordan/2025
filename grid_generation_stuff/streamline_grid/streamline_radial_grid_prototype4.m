%% streamline/radial_grid prototype 4 (07/13/2026)
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

airfoil_inputs = struct();
airfoil_inputs.epsilon = 0.1;
airfoil_inputs.kappa   = 0.1;
airfoil_inputs.tau     = 0.0;
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = 0; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;
boundary_distance = 500;
n_body_pts = 129;
n_wake_pts = 65;
n_transition = 1;
jmax_in    = 1;
AR_LE      =  1.0;
AR_TE      =  1.0;
delta_LE   = 0.001;
airfoil = local_airfoil_generator(airfoil_inputs);
[x,y,theta0,F1,psi_1,phi_TE,h_LE,h_TE] = get_airfoil_streamline_body_coordinates( airfoil, ...
                                                         n_body_pts, ...
                                                         n_transition, ...
                                                         delta_LE,     ...
                                                         AR_LE, AR_TE );

dj = 50;
di = 65;

hold on
% plot(x,y,'r.-')
axis equal

zeta1 = airfoil.zeta_from_theta(theta0);
plot(real(zeta1),imag(zeta1),'k')

z1 = airfoil.airfoil_coords(theta0);
plot(real(z1),imag(z1),'k');

zeta2 = airfoil.zeta_from_z(z1);
plot(real(zeta2),imag(zeta2),'r--');

zeta2 = airfoil.zeta_from_z_alt(z1);
plot(real(zeta2),imag(zeta2),'b--');

z2 = airfoil.z_from_zeta(zeta2);
plot(real(z2),imag(z2),'r--');

zeta3 = airfoil.zeta_from_theta_r(theta0,1.2*airfoil.a);
plot(real(zeta3),imag(zeta3),'b')

z3 = airfoil.z_from_zeta(zeta3);
plot(real(z3),imag(z3),'b');

zeta4 = airfoil.zeta_from_z(z3);
plot(real(zeta4),imag(zeta4),'g--');

z4 = airfoil.z_from_zeta(zeta4);
plot(real(z4),imag(z4),'g--');

function airfoil = local_airfoil_generator(in)
% epsilon = 0.1;
% kappa   = 0.0;
% tau     = 0.0;
% vinf    = 75.0;
% rhoinf  = 1.0;
% pinf    = 100000.0;
% gamma   = 1.4;
% alpha   = 0; % (degrees)
% rho_ref = 1.0;
% p_ref   = 100000.0;
a_ref   = sqrt(in.gamma*in.p_ref/in.rho_ref);

% nondimensionalize inputs
vinf   = in.vinf/a_ref;
rhoinf = in.rhoinf/in.rho_ref;
pinf   = in.pinf/(in.rho_ref*a_ref^2);

airfoil        = kt_airfoil( in.epsilon, in.kappa, in.tau );
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(in.alpha);
end