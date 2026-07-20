%% looking at airfoil potential field (07/19/2026)
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
airfoil_inputs.alpha   = 5; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;


boundary_distance = 0.5;
n_body_pts        = 257;
n_transition      = 1;
jmax              = 101;

AR_LE             = 1.0;
AR_TE             = 1.0;
delta_LE          = 0.001;
airfoil = local_airfoil_generator(airfoil_inputs);
[ theta_u,theta_l, psi_1, psi_FF, phi_TE, h_LE, h_TE, delta_TE ] = get_streamline_geometry( ...
                                                            airfoil,  ...
                                                            delta_LE, ...
                                                            AR_LE,    ...
                                                            AR_TE,    ...
                                                            boundary_distance );
[x,y,theta,phi,f1,df1,ddf1] = get_airfoil_coordinates( airfoil, delta_LE, delta_TE, n_body_pts, n_transition );
r_vec = my_geomspace(jmax,0,xmax=boundary_distance,dx0=h_LE);
r_vec = airfoil.a + airfoil.chord*r_vec;
[THETA,R] = ndgrid(theta,r_vec);
ZETA1 = airfoil.zeta_from_theta_r(THETA,R);
Z1 = airfoil.zs_from_zeta(ZETA1);
X1 = real(Z1); Y1 = imag(Z1);
F0 = airfoil.F_cylinder(ZETA1);
contourf(X1,Y1,real(F1))
% contourf(real(ZETA1),imag(ZETA1),imag(airfoil.w_cylinder(ZETA1)))
axis equal
% view(2)
colorbar

hold on
zsp = airfoil.zs_from_theta( airfoil.thetaSP);
plot(real(zsp),imag(zsp),'rx')

% csape([] {'complete','complete'}


function [theta_u,theta_l,psi_1,psi_FF,phi_TE,h_LE,h_TE,delta_TE] = get_streamline_geometry( airfoil, delta_LE, AR_LE, AR_TE, boundary_distance )
% given streamwise spacing at LE, target AR (streamwise/normal)
% for LE and TE, need to find trailing edge spacing

% approximate cell height at LE
h_LE = delta_LE/AR_LE;

% compute the angles where phi=0

% upper surface
theta_u = fzero( @(theta) real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta) ) ), [0,pi] );
z_u = airfoil.z_from_zeta(airfoil.zeta_from_theta(theta_u));
z_u = (z_u-airfoil.zLE)/airfoil.chord;
z_u = z_u + h_LE*airfoil.unit_normal_cmplx(theta_u);
% find the corresponding value of psi
psi_u = airfoil.psi_from_xy( real(z_u), imag(z_u) );

% lower surface
theta_l = fzero( @(theta) real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta) ) ), [pi,2*pi] );
z_l = airfoil.z_from_zeta(airfoil.zeta_from_theta(theta_l));
z_l = (z_l-airfoil.zLE)/airfoil.chord;
z_l = z_l + h_LE*airfoil.unit_normal_cmplx(theta_l);
% find the corresponding value of psi
psi_l = airfoil.psi_from_xy( real(z_l), imag(z_l) );

% let's average them to maybe get the height at LE closer to target ...
psi_1 = sign(psi_u)*0.5*( abs(psi_u) + abs(psi_l) );

% from here to the leading edge, we are assuming approximately constant cell height

% compute the height at the trailing edge for this psi
zTE = airfoil.zTE;
zTE = (zTE-airfoil.zLE)/airfoil.chord;
phi_TE = airfoil.phi_from_xy(real(zTE),imag(zTE));
z3 = airfoil.z_from_phi_psi(phi_TE,psi_1);
h_TE = abs(z3-zTE);

% compute the corresponding surface spacing
delta_TE = AR_TE*h_TE;

% compute outer boundary streamline

% construct a circle, centered at the trailing edge, and find a location
% such that it has the same potential (phi) as at the trailing edge
%       The termination of this potential line will then be exactly
%       boundary distance away from the trailing edge
%       (at least at that terminal point)

xfun = @(theta) real(zTE) + boundary_distance.*cos(theta);
yfun = @(theta) imag(zTE) + boundary_distance.*sin(theta);
objfun = @(theta) phi_TE-airfoil.phi_from_xy(xfun(theta),yfun(theta));
theta = fzero(@(theta)objfun(theta),pi/2);
x_ff = xfun(theta);
y_ff = yfun(theta);
psi_FF = airfoil.psi_from_xy(x_ff,y_ff);
end

function [x,y,theta,phi,F1,dF1,ddF1] = get_airfoil_coordinates( airfoil, delta_LE, delta_TE, n_pts, n_transition )
% get body coordinates, angle in zeta plane, and stretching function

% get parametric location of maximum curvature on the airfoil
% which will approximately be the leading edge
% tLE = airfoil.thetaCmax/(2*pi);

tLE = airfoil.thetaLE/(2*pi);

% arc length
L = airfoil.airfoil_arc_length(0,1);

% get normalized arc length to point of max curvature
t0 = airfoil.airfoil_arc_length(0,tLE)/L;

% nondimensionalize spacings
cxL   = airfoil.chord/L;

dLE = delta_LE*cxL;
dTE = delta_TE*cxL;

% get buffer distance over which to blend the two functions
offset = ( (n_transition-1)/(n_pts-1) );

% generate asymmetric stretching function for airfoil surface
[F1,dF1,ddF1] = hermite_blend_2_vinokur_asym(n_pts,t0,dLE,dTE,offset,true);
[x,y,theta] = airfoil.output_airfoil_coords1(n_pts,F1);

% flip to match the clockwise convention
x = flip(x);
y = flip(y);
theta = flip(theta);

% potential along airfoil
phi   = real(airfoil.F_cylinder(airfoil.zeta_from_theta(theta)));
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