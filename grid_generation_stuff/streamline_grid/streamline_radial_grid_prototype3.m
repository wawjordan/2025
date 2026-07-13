%% streamline/radial_grid prototype 3 (07/12/2026)
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
airfoil_inputs.tau     = 0.1;
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

dj = 10;
di = 65;
overlap = 0;

hold on
plot(x,y,'r.-')
axis equal

% find phi at trailing edge, then find coordinates for that phi at the outer boundary
[phi1,x1,y1,~] = min_phi_on_circle(airfoil,x(1),y(1),boundary_distance,pi/2);

% get the streamline that goes through that point
psi_jmax = airfoil.psi_from_xy(x1,y1);

% phi0 = airfoil.phi_from_xy(x,y);
phi0 = airfoil.phi_from_xy(x,y);

phi02 = airfoil.F_cylinder(airfoil.zeta_from_theta(theta0));
% phi = phi0(1:di);


X = zeros(di,dj);
Y = zeros(di,dj);

[~,idx] = min(abs(phi0));
idx = 20;
idx = idx+1;
phi = phi0(1:idx+overlap);
psi0 = my_geomspace(129,0,xmax=psi_jmax,dx0=psi_1);
psi = psi0(1:dj);
[PHI,PSI] = ndgrid(phi,psi);
Z1 = airfoil.z_from_phi_psi(PHI,PSI);
X1 = real(Z1); Y1 = imag(Z1);
plot(X1,Y1,'k')
plot(X1.',Y1.','k')
X(1:idx,:) = X1;
Y(1:idx,:) = Y1;



z1 = airfoil.z_from_phi_psi(phi(idx-overlap),psi);
plot(real(z1),imag(z1),'g')
[theta1,r] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z1+airfoil.zLE));
[~,idx1] = min(abs(theta0-airfoil.thetaCmax));
[~,idx2] = min(abs(theta0-theta1(1)));


theta = theta0(idx1:idx2);
[THETA,R] =ndgrid(theta,r);
Z2 = airfoil.z_from_theta_r(THETA,R);
X2 = real(Z2); Y2 = imag(Z2);
plot(X2,Y2,'b')
plot(X2.',Y2.','b')

function [phi,x,y,theta] = min_phi_on_circle(airfoil,x0,y0,r,theta_guess)
xfun = @(theta) x0 + r.*cos(theta);
yfun = @(theta) y0 + r.*sin(theta);
phi0 = airfoil.phi_from_xy(x0,y0);
objfun = @(theta) phi0-airfoil.phi_from_xy(xfun(theta),yfun(theta));
theta = fzero(@(theta)objfun(theta),theta_guess);
x = xfun(theta);
y = yfun(theta);
phi = airfoil.phi_from_xy(x,y);
end

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