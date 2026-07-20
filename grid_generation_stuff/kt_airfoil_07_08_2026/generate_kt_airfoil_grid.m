function [OUT,airfoil] = generate_kt_airfoil_grid(airfoil_inputs,n_body_pts,n_wake_pts,jmax,AR_LE,AR_TE,delta_LE,delta_TE,spline_order)
% number of points off body
% number of points on the body
% number of points in wake
% target aspect ratio at the leading edge
% target aspect ratio at the trailing edge
% target wall spacing at leading edge
% target wall spacing at trailing edge

OUT = struct();
airfoil = local_airfoil_generator(airfoil_inputs);
n_transition = 1;   % number of transition points for surface spacing blending
imax = n_body_pts + 2*(n_wake_pts-1);
i1 = n_wake_pts;
i2 = n_wake_pts+n_body_pts-1;
wake_multiplier = 1;
boundary_distance = 500;
mu                = 0.1;             % 4th order explicit smoothing factor
muim              = 0.5;             % Implicit smoothing factor
scjmax = bivariate_ramps(imax,jmax,[1.0000,0.999,1.0,1.0,0.999,1.0000],[0,1/3,15/32,17/32,2/3,1],...
                                   [1.0000,1.0000,0.98],[0,0.1,1]);
alpham = bivariate_ramps(imax,jmax,[1.0000,2.0000,2.0000,1.0000],[0,3/8,5/8,1],...
                                   [0.5000,1.0000,2.0000,0.5000],[0,1/8,3/8,1]);
[x_airfoil,y_airfoil,theta,theta_fun,~,~,~,~,fy_wake] = get_airfoil_coordinates(airfoil,n_body_pts,n_transition,AR_LE,AR_TE,delta_LE,delta_TE,n_wake_pts,boundary_distance);
[x,y] = hyperbolic_C_grid_local_v2( x_airfoil, y_airfoil,                 ...
                                      n_wake_pts, jmax, ...
                                      boundary_distance=boundary_distance, ...
                                      wake_multiplier=wake_multiplier, ...
                                      AR_LE = AR_LE, ...
                                      AR_TE = AR_TE, ...
                                      n_extrude_layers=0, ...
                                      n_radial_passes=1,...
                                      scjmax = scjmax, ...
                                      mu = mu, ...
                                      muim = muim, ...
                                      alpham = alpham, ...
                                      fy_wake = fy_wake );

% figure()
% s = 2^0;
% hold on
% plot(x(1:s:end,1:s:end),y(1:s:end,1:s:end),'k');
% plot(x(1:s:end,1:s:end).',y(1:s:end,1:s:end).','k');
% axis equal
OUT.airfoil_inputs = airfoil_inputs;
OUT.base_grid = struct();
OUT.base_grid.n_body_pts = n_body_pts;
OUT.base_grid.n_wake_pts = n_wake_pts;
OUT.base_grid.jmax       = jmax;
OUT.base_grid.imax       = imax;
OUT.base_grid.i1         = i1;
OUT.base_grid.i2         = i2;
OUT.base_grid.AR_LE      = AR_LE;
OUT.base_grid.AR_TE      = AR_TE;
OUT.base_grid.delta_LE   = delta_LE;
OUT.base_grid.delta_TE   = delta_TE;
OUT.spline_order         = spline_order;
OUT.base_grid.boundary_distance = boundary_distance;
OUT.base_grid.x_airfoil  = x_airfoil;
OUT.base_grid.y_airfoil  = y_airfoil;
OUT.base_grid.x          = x;
OUT.base_grid.y          = y;
OUT.base_grid.theta      = theta;
OUT.base_grid.theta_fun  = theta_fun;
if (jmax>1)
OUT.xs = spapi({spline_order,spline_order},{linspace(0,1,imax),linspace(0,1,jmax)},x);
OUT.ys = spapi({spline_order,spline_order},{linspace(0,1,imax),linspace(0,1,jmax)},y);
end
OUT.xi_min  = 0;
OUT.xi_max  = 1;
OUT.eta_min = 0;
OUT.eta_max = 1;

end

function f = bivariate_ramps(imax,jmax,i_vals,i_breaks,j_vals,j_breaks)
f_i = @(i) smooth_ramps((i-1)/(imax-1),i_breaks,i_vals);
f_j = @(j) smooth_ramps((j-1)/(jmax-1),j_breaks,j_vals);
f = @(i,j) f_i(i).*f_j(j);
end

function [x,y,theta,theta_fun,F1,dF1,ddF1,fx_wake,fy_wake] = get_airfoil_coordinates( ...
                                             airfoil, ...
                                             n_body_pts, ...
                                             n_transition, ...
                                             AR_LE,        ...
                                             AR_TE,        ...
                                             delta_LE,     ...
                                             delta_TE,     ...
                                             n_wake_pts,   ...
                                             boundary_distance )
% get parametric location of maximum curvature on the airfoil
% which will approximately be the leading edge
tLE = airfoil.thetaCmax/(2*pi);
% tLE = airfoil.thetaSP/(2*pi);


L = airfoil.airfoil_arc_length(0,1);

% get normalized arc length to point of max curvature
t0 = airfoil.airfoil_arc_length(0,tLE)/L;

% nondimensionalize spacings
% cxL   = airfoil.chord/L;
% perim = L/airfoil.chord;
% uniform_spacing = perim/(n_body_pts-1);
% dLE = AR_LE*delta_LE*uniform_spacing*cxL;
% dTE = AR_TE*delta_TE*uniform_spacing*cxL;

uniform_spacing = 1/(n_body_pts-1);
dLE = AR_LE*delta_LE*uniform_spacing;
dTE = AR_TE*delta_TE*uniform_spacing;

% dLE = delta_LE*cxL;
% dTE = delta_TE*cxL;

% get buffer distance over which to blend the two functions
% offset = ( (n_transition-1)/(n_body_pts-1) )*cxL;
offset = ( (n_transition-1)/(n_body_pts-1) );

% generate asymmetric stretching function for airfoil surface
[F1,dF1,ddF1] = hermite_blend_2_vinokur_asym(n_body_pts,t0,dLE,dTE,offset,true);
% % [f,df,ddf] = hermite_blend_2_vinokur_asym(N,t0,d0,d1,off,refine)
[x,y,theta] = airfoil.output_airfoil_coords1(n_body_pts,F1);

f1 = @(theta) abs(airfoil.airfoil_differential_arc_length(theta));
f1c = chebfun(f1,[0,2*pi],'splitting','on');
s_c = cumsum(f1c)./sum(f1c,0,2*pi);
theta_s = inv(s_c,'splitting','on');
theta_fun = @(t) theta_s( F1(1-t) );
% flip to match the clockwise convention
x = flip(x);
y = flip(y);
theta = flip(theta);
if (nargin > 8 && n_wake_pts>0)
    [ ~, ~, fx_wake ] = wake_cut_pts( x, y, boundary_distance, n_wake_pts );
    fy_wake = @(t) airfoil.psi_y_from_x(0,fx_wake(t),0);
else
    fx_wake = [];
    fy_wake = [];
end

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