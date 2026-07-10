%% airfoil grid gen prototyping (05/11/2026)
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
% y_from_x = chebfun(@(x) airfoil.streamline_y(0,x,0),[1,500],1000);
% s_from_x = cumsum(sqrt(1+diff(y_from_x).^2));
% x_from_s = inv(s_from_x);
% L = s_from_x(500);
% y_from_t = @(t) y_from_x( x_from_s(t*L) );
% y_from_t = @(t) y_from_x( 1 + (500-1)*t);
% pinf = minimax(f, deg, 'tol', 1e-16);
% p2   = polyfit(f,deg);
% r1   = aaa(f);

% perim = airfoil.airfoil_arc_length(0,1)/airfoil.chord;
% uniform_spacing = perim/(n_body_pts-1);
% delta_LE = 0.0005; % target wall spacing at leading edge
% delta_TE = 0.0005; % target wall spacing at trailing edge
% delta_LE = 0.1*uniform_spacing; % target wall spacing at leading edge
% delta_TE = 0.1*uniform_spacing; % target wall spacing at trailing edge
% [x_airfoil,y_airfoil,F1,dF1,ddF1] = get_airfoil_coordinates(airfoil,n_body_pts,n_transition,delta_LE*AR_LE,delta_TE*AR_TE);
% [x2,y2] = hyperbolic_C_grid_local_v2( x_airfoil, y_airfoil,                 ...
%                                  boundary_distance, jmax, ...
%                                  AR_LE, AR_TE, scjmax, ...
%                                  mu, muim, alpham, ...
%                                  n_wake_pts, wake_multiplier );

% jmax         = 513;  % number of points off body
% n_body_pts   = 1025;  % number of points on the body
% n_wake_pts   = 257;  % number of points in wake
% n_transition = 33;   % number of transition points for surface spacing blending
jmax         = 129;  % number of points off body
n_body_pts   = 129;  % number of points on the body
n_wake_pts   = 129;  % number of points in wake
n_transition = 1;   % number of transition points for surface spacing blending
imax = n_body_pts + 2*(n_wake_pts-1);
i1 = n_wake_pts;
i2 = n_wake_pts+n_body_pts-1;

AR_LE    = 1; % target aspect ratio at the leading edge
AR_TE    = 1; % target aspect ratio at the trailing edge
delta_LE = 0.1; % target wall spacing at leading edge
delta_TE = 0.1; % target wall spacing at trailing edge
wake_multiplier = 1;
boundary_distance = 500;
mu                = 0.1;             % 4th order explicit smoothing factor
muim              = 0.5;             % Implicit smoothing factor
scjmax            = 0.997;           % Scaling Factor
alpham            = 4;               % Alpha scheme integration factor
% o = 33;
% f1s =  @(t) smooth_transition_1_side(t,i2-o,imax,0,1-0.9999);
% f2s =  @(t) smooth_transition_1_side(t,1,i1+o,1,0.9999);
% f1ss = @(t) smooth_transition_1_side(t,1,jmax,1,0.99);
% scjmax = @(i,j) (f1s(i) + f2s(i)).*f1ss(j);
% f1a = @(t) smooth_transition_1_side(t,1,jmax/3,0.5,3);
% f2a = @(t) smooth_transition_1_side(t,  jmax/3,2*jmax/3,0,0.5-3);
% alpham = @(i,j) f1a(j) + f2a(j) + 0*i;

% scjmax = bivariate_ramps(imax,jmax,[1.0000,0.999,0.999,1.0000],[0,1/3,2/3,1],...
%                                    [1.0000,0.9950],[0,1]);
% alpham = bivariate_ramps(imax,jmax,[1.0000,2.0000,2.0000,1.0000],[0,3/8,5/8,1],...
%                                    [0.5000,1.0000,2.0000,0.5000],[0,1/8,3/8,1]);

scjmax = bivariate_ramps(imax,jmax,[1.0000,0.999,1.0,1.0,0.999,1.0000],[0,1/3,15/32,17/32,2/3,1],...
                                   [1.0000,1.0000,0.98],[0,0.1,1]);
alpham = bivariate_ramps(imax,jmax,[1.0000,2.0000,2.0000,1.0000],[0,3/8,5/8,1],...
                                   [0.5000,1.0000,2.0000,0.5000],[0,1/8,3/8,1]);



[x_airfoil,y_airfoil,F1,dF1,~,~,fy_wake] = get_airfoil_coordinates(airfoil,n_body_pts,n_transition,AR_LE,AR_TE,delta_LE,delta_TE,n_wake_pts,boundary_distance);


fplot(dF1,[0.45,0.55])

hold on
plot(x_airfoil,y_airfoil,'r.-');
axis equal

x1 = (airfoil.xTE-airfoil.xLE)/airfoil.chord;
y1 = airfoil.yTE/airfoil.chord;
[xfun,yfun] = potential_grid_line_v1(airfoil,x1,y1);
psi_bnd = fzero(@(psi)yfun(psi)-boundary_distance,4*boundary_distance);
tpsi = linspace(-1,1,101);
plot(xfun(tpsi),yfun(tpsi));
% find phi for Trailing edge
% phiTE0 = airfoil.phi_from_xy( x1,y1);
% ty = linspace(-500,500,10001);
% plot(airfoil.phi_x_from_y(phiTE0,ty,0),ty,'b')
% xlim([-1,2])

tx = linspace(0,500,10001);
% plot(tx,airfoil.psi_y_from_x(500*airfoil.chord,tx,-0.1),'g')
% plot(airfoil.phi_x_from_y(500*airfoil.chord,ty,0),ty,'b')


tic;
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
time = toc;
figure()
s = 2^0;
hold on
plot(x(1:s:end,1:s:end),y(1:s:end,1:s:end),'k');
plot(x(1:s:end,1:s:end).',y(1:s:end,1:s:end).','k');
axis equal


spline_x = spapi({5,5},{linspace(0,1,imax),linspace(0,1,jmax)},x);
spline_y = spapi({5,5},{linspace(0,1,imax),linspace(0,1,jmax)},y);


ni = 1/4;
nj = 1/4;
% ni = 8;
% nj = 8;
x2 = fnval(spline_x,{linspace(0,1,(imax-1)*ni+1),linspace(0,1,(jmax-1)*nj+1)});
y2 = fnval(spline_y,{linspace(0,1,(imax-1)*ni+1),linspace(0,1,(jmax-1)*nj+1)});

figure()
plot_edge_length_ratio(x2(1:s:end,1:s:end),y2(1:s:end,1:s:end),1,'jet',1);
xlim([-1,2])
ylim([-1.5,1.5])

plot_edge_length_ratio(x2(1:s:end,1:s:end),y2(1:s:end,1:s:end),2,'jet',1);
xlim([-1,2])
ylim([-1.5,1.5])

% L1 = edge_length_ratio_2D(x2,y2,1);
% L2 = edge_length_ratio_2D(x2,y2,2);

hold off


function hfig = plot_edge_length_ratio(x,y,dir,cmap,linewidth)
hfig = figure();
set(gca,'Color','k');
colormap(cmap);
L = edge_length_ratio_2D(x,y,dir);
x = padarray(x,[1,1],nan,'both');
y = padarray(y,[1,1],nan,'both');
L = padarray(L,[1,1],nan,'both');
hold on;
patch(x,y,L,'edgecolor','interp','Marker','.','LineWidth',linewidth);
patch(x.',y.',L.','edgecolor','interp','LineWidth',linewidth);
axis equal
colorbar;
hold off
end

function f = bivariate_ramps(imax,jmax,i_vals,i_breaks,j_vals,j_breaks)
f_i = @(i) smooth_ramps((i-1)/(imax-1),i_breaks,i_vals);
f_j = @(j) smooth_ramps((j-1)/(jmax-1),j_breaks,j_vals);
f = @(i,j) f_i(i).*f_j(j);

end

function [x,y,F1,dF1,ddF1,fx_wake,fy_wake] = get_airfoil_coordinates( ...
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
[x,y] = airfoil.output_airfoil_coords1(n_body_pts,F1);

% flip to match the clockwise convention
x = flip(x);
y = flip(y);
if (nargin > 8 && n_wake_pts>0)
    [ ~, ~, fx_wake ] = wake_cut_pts( x, y, boundary_distance, n_wake_pts );
    fy_wake = @(t) airfoil.psi_y_from_x(0,fx_wake(t),0);
else
    fx_wake = [];
    fy_wake = [];
end

end

function airfoil = local_airfoil_generator()
epsilon = 0.1;
kappa   = 0.0;
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