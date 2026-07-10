%% streamline/radial_grid (07/09/2026)
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
AR_LE      =  10.0;
AR_TE      =  10.0;
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
[phi1,x1,y1,theta] = min_phi_on_circle(airfoil,x0,y0,r,theta_guess);

x_front = @(y) airfoil.phi_x_from_y(phi1,y,x0);

% get the streamline that goes through that point
psi1 = airfoil.psi_from_xy(x1,y1);

y_upper = @(x) airfoil.psi_y_from_x(psi1,x,y1);

% find where this intercepts the potential line eminating from the last wake pt
x2 = x(1);
y2 = y(1);
phi2 = airfoil.phi_from_xy(x2,y2);
x_back  = @(y) airfoil.phi_x_from_y(phi2,y,x2);
psifun = @(y) airfoil.psi_from_xy(x_back(y),y);
options = optimset('TolFun',1e-15,'TolX',1e-17);
obj_fun = @(y) psifun(y) - psi1;
y3 = fzero(@(y)obj_fun(y),y_upper(x2),options);
x3 = x_back(y3);

plot(x0,y0,'ro');
plot(x1,y1,'bo');
plot(x2,y2,'go');
plot(x3,y3,'mo');

y_front = linspace(y0,y1,101);
plot(x_front(y_front),y_front)

x_upper = linspace(x1,x3,101);
plot(x_upper,y_upper(x_upper))

y_back = linspace(y2,y3,101);
plot(x_back(y_back),y_back)

% then iterate over the potential lines for the upper surface of the airfoil


% imax = OUT.base_grid.imax;
i2 = OUT.base_grid.i2;
i1 = i2 - 20;
for i = i1:i2
    % find where this intercepts the potential line eminating from the last wake pt
    xi = x(i);
    yi = y(i);
    phi_i = airfoil.phi_from_xy(xi,yi);
    x_fi  = @(y) airfoil.phi_x_from_y(phi_i,y,xi);
    psifun = @(y) airfoil.psi_from_xy(x_fi(y),y);
    obj_fun = @(y) psifun(y) - psi1;
    yi_ = fzero(@(y)obj_fun(y),y_upper(xi),options);
    % xi_ = x_fi(y3);
    y_tmp = linspace(yi,yi+0.5,101);
    plot(x_fi(y_tmp),y_tmp)
end
z0 = 1i;
% phi = airfoil.phi_from_xy()
z1 = airfoil.z_from_phi_psi(linspace(0,5,101),0.1);
plot(real(z1),imag(z1),'g')


% function [x,y] = airfoil_radial_potential_grid(airfoil,x_airfoil,y_airfoil,boundary_distance)
% % find potential at location (x,y)
% phi  = airfoil.phi_from_xy(x,y);
% psi0 = airfoil.psi_from_xy(x,y);
% xfun1  = @(y) airfoil.phi_x_from_y(phi,y,x);
% psifun = @(y) airfoil.psi_from_xy(xfun1(y),y);
% obj_fun = @(psi,y) psifun(y) - psi;
% yfun = @(psi) arrayfun(@(psi)fzero(@(y)obj_fun(psi,y),psi0,options),psi);
% xfun = @(psi) xfun1(yfun(psi));
% end
% 
% function [x,y] = wake_pts(airfoil,xjmin,yjmin,boundary_distance)
% 
% 
% end

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

function z = stream_coordinate(airfoil,phi,psi,z0,max_iter)

airfoil.F_cylinder( airfoil.zeta_from_z)
res0 = abs()




end