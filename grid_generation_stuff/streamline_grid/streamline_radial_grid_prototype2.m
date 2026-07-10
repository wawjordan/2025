%% streamline/radial_grid prototype 2 (07/09/2026)
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
airfoil_inputs.kappa   = 0.0;
airfoil_inputs.tau     = 0.0;
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = 0; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;

n_body_pts = 129;
n_wake_pts = 65;
jmax_in    = 1;
AR_LE      =  1.0;
AR_TE      =  1.0;
delta_LE   = 0.1;
delta_TE   = 0.1;
spline_order = 5;
[OUT,airfoil] = generate_kt_airfoil_grid(airfoil_inputs,n_body_pts,n_wake_pts,jmax_in,AR_LE,AR_TE,delta_LE,delta_TE,spline_order);
x = OUT.base_grid.x(:,1);
y = OUT.base_grid.y(:,1);
hold on
plot(x,y,'k')
axis equal

% find phi at trailing edge, then find coordinates for that phi at the outer boundary
theta_guess = pi/2;
r  = OUT.base_grid.boundary_distance;
x0 = OUT.base_grid.x_airfoil(1);
y0 = OUT.base_grid.y_airfoil(1);
[phi1,x1,y1,~] = min_phi_on_circle(airfoil,x0,y0,r,theta_guess);

% get the streamline that goes through that point
psi1 = airfoil.psi_from_xy(x1,y1);

phi0 = airfoil.phi_from_xy(x,y);
mask = phi0>=0;
phi = phi0(mask);
psi = my_geomspace(129,0,xmax=psi1,dx0=0.1);
[PHI,PSI] = ndgrid(phi,psi);
Z = airfoil.z_from_phi_psi(PHI,PSI);
X = real(Z); Y = imag(Z);
plot(X,Y,'k')
plot(X.',Y.','k')


z1 = airfoil.z_from_phi_psi(0,psi);
[theta0,r] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z1+airfoil.zLE));

z2 = airfoil.z_from_theta_r(theta0,r);
plot(real(z1),imag(z1),'k')
plot(real(z2),imag(z2),'r--')


theta = OUT.base_grid.theta;
mask = theta>=pi/2 & theta<=3*pi/2;


[R,THETA] =ndgrid(r,theta(mask));
Z2 = airfoil.z_from_theta_r(THETA,R);
X2 = real(Z2); Y2 = imag(Z2);
plot(X2,Y2,'k')
plot(X2.',Y2.','k')


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