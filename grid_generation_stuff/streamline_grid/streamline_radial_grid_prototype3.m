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
airfoil_inputs.kappa   = 0.0;
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

dj = 30;
di = n_body_pts;

hold on
plot(x,y,'r.-')
axis equal

% find phi at trailing edge, then find coordinates for that phi at the outer boundary
[phi1,x1,y1,~] = min_phi_on_circle(airfoil,x(1),y(1),boundary_distance,pi/2);
i = 65;
% x0 =  x(end);
% y0 =  y(end)-1e-2;
% [xs,ys,ts] = ode_streamline(airfoil,x0,y0,2000,false);
% plot(xs,ys,'b')
% [xs,ys,ts] = ode_streamline(airfoil,x0,y0,2000,true);
% plot(xs,ys,'b')

% get the streamline that goes through that point
psi_jmax = airfoil.psi_from_xy(x1,y1);

psi0 = my_geomspace(129,0,xmax=psi_jmax,dx0=psi_1);
% phi0 = airfoil.phi_from_xy(x,y);
phi0 = real(airfoil.F_cylinder(airfoil.zeta_from_theta(theta0)));
X = zeros(di,dj);
Y = zeros(di,dj);

X(:,1) = x;
Y(:,1) = y;


[~,idx1] = min(abs(theta0-(airfoil.thetaSP+pi/2)));
[~,idx2] = min(abs(theta0-(airfoil.thetaSP-pi/2)));

%% lower surface
phi =  phi0(1:idx1);
psi = -psi0(2:dj);
[PHI,PSI] = ndgrid(phi,psi);
Z1 = airfoil.z_from_phi_psi(PHI,PSI);
X1 = real(Z1); Y1 = imag(Z1);
X(1:idx1,2:dj) = X1; Y(1:idx1,2:dj) = Y1;
% plot(X1,Y1,'k')
% plot(X1.',Y1.','k')


%% leading edge
z1 = airfoil.z_from_phi_psi(phi0(idx1),psi);
plot(real(z1),imag(z1),'g')
[theta1,r1] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z1+airfoil.zLE));
phi =  phi0(idx2:end);
psi =  psi0(2:dj);
z2 = airfoil.z_from_phi_psi(phi0(idx2),psi);
plot(real(z2),imag(z2),'g')
[theta2,r2] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z2+airfoil.zLE));


theta = theta0(idx1:idx2);
[THETA,R] =ndgrid(theta,r1);
for j = 1:numel(r1)
    R(:,j) = linspace(r1(j),r2(j),numel(theta));
    % for i = 1:numel(theta)
    %     R(i,j) = airfoil.r_from_psi_theta(psi(j),theta(i));
    % end
end
Z2 = airfoil.z_from_theta_r(THETA,R);
X2 = real(Z2); Y2 = imag(Z2);
X(idx1:idx2,2:dj) = X2; Y(idx1:idx2,2:dj) = Y2;
% plot(X2,Y2,'b')
% plot(X2.',Y2.','b')

%% upper surface
[PHI,PSI] = ndgrid(phi,psi);
Z3 = airfoil.z_from_phi_psi(PHI,PSI);
X3 = real(Z3); Y3 = imag(Z3);
X(idx2+1:end,2:dj) = X3(2:end,:); Y(idx2+1:end,2:dj) = Y3(2:end,:);
% plot(X3,Y3,'k')
% plot(X3.',Y3.','k')

plot(X,Y,'k')
plot(X.',Y.','k')

idx1;

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

function [x,y,ts] = ode_streamline(airfoil,x0,y0,t,potential)
tspan = [0 t];
z0 = x0 + 1i*y0;
z0 = z0*airfoil.chord + airfoil.zLE;
zeta0 = airfoil.zeta_from_z_alt(z0);
ys0 = [real(zeta0), imag(zeta0)];
w0  = airfoil.w_cylinder(zeta0);
options = odeset('RelTol',1e-12,'AbsTol',1e-16,'InitialStep',1e-12,'NormControl','on');
if ( potential )
    yp0 = [imag(w0),real(w0)];
    options = odeset(options,'InitialSlope',yp0);
    % phi0 = airfoil.phi_from_xy(x0,y0);
    % s = 2*(phi0>=0)-1;
    psi0 = airfoil.psi_from_xy(x0,y0);
    s = 2*(psi0>=0)-1;
    [ts,ys] = ode113(@(t,y) phi_integrate(t,y,airfoil), s*tspan, ys0, options);
    ts = s*ts;
else
    yp0 = [real(w0),-imag(w0)];
    options = odeset(options,'InitialSlope',yp0);
    [ts,ys] = ode113(@(t,y) psi_integrate(t,y,airfoil), tspan, ys0, options);
end
zeta = ys(:,1) + 1i*ys(:,2);
z = airfoil.z_from_zeta(zeta);
z = (z-airfoil.zLE)/airfoil.chord;
x = real(z);
y = imag(z);
end

function dydt = phi_integrate(~,y,airfoil)
    dydt = psi_integrate(0,y,airfoil);
    tmp = dydt(1);
    dydt(1) = -dydt(2);
    dydt(2) = tmp;
end

function dydt = psi_integrate(~,y,airfoil)
zeta = y(1) + 1i*y(2);
% w = airfoil.w_airfoil(zeta);
w = airfoil.w_cylinder(zeta);
dydt = zeros(2,1);
dydt(1) = real(w);
dydt(2) =-imag(w);
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