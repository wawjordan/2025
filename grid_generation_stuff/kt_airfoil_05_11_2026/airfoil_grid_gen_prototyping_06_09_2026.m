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
% n_body_pts   = 513;  % number of points on the body
% n_wake_pts   = 257;  % number of points in wake
% n_transition = 33;   % number of transition points for surface spacing blending
jmax         = 65;  % number of points off body
n_body_pts   = 65;  % number of points on the body
n_wake_pts   = 33;  % number of points in wake
n_transition = 9;   % number of transition points for surface spacing blending
imax = n_body_pts + 2*(n_wake_pts-1);
i1 = n_wake_pts;
i2 = n_wake_pts+n_body_pts-1;

AR_LE    = 1; % target aspect ratio at the leading edge
AR_TE    = 1; % target aspect ratio at the trailing edge
delta_LE = 0.05; % target wall spacing at leading edge
delta_TE = 0.05; % target wall spacing at trailing edge
wake_multiplier = 1;
boundary_distance = 500;
% scjmax            = 0.99;           % Scaling Factor
mu                = 0.1;             % 4th order explicit smoothing factor
muim              = 0.5;             % Implicit smoothing factor
% scjmax            = 0.9;           % Scaling Factor
o = 33;
% f1s = @(t) smooth_transition_1_side(t,i2,imax+1-o,0,1-0.999);
% f2s = @(t) smooth_transition_1_side(t,o,i1,1,0.999);
% f1ss = @(t) smooth_transition_1_side(t,1,jmax,1,0.99);
f1s = @(t) smooth_transition_1_side(t,i2-o,imax,0,1-0.9999);
f2s = @(t) smooth_transition_1_side(t,1,i1+o,1,0.9999);
f1ss = @(t) smooth_transition_1_side(t,1,jmax,1,0.99);
scjmax = @(i,j) (f1s(i) + f2s(i)).*f1ss(j);
% f1s = @(t) smooth_transition_1_side(t,1,jmax,1,0.99);
% scjmax = @(i,j) f1s(j) + 0*i;
% alpham            = 0.5;            % Alpha scheme integration factor
% f1a = @(t) smooth_transition_1_side(t,1,jmax/3,0.5,3);
% f2a = @(t) smooth_transition_1_side(t,  jmax/3,2*jmax/3,0,0.5-1);
f1a = @(t) smooth_transition_1_side(t,1,jmax/6,0.5,3);
f2a = @(t) smooth_transition_1_side(t,  jmax/6,2*jmax/3,0,0.5-1);
% f2a =  @(t) 0*t;
alpham = @(i,j) f1a(j) + f2a(j) + 0*i;




[x_airfoil,y_airfoil,~,~,~,~,fy_wake] = get_airfoil_coordinates(airfoil,n_body_pts,n_transition,AR_LE,AR_TE,delta_LE,delta_TE,n_wake_pts,boundary_distance);




hold on
plot(x_airfoil,y_airfoil,'r.-');
axis equal

% phiTE0 = airfoil.phi_from_xy(airfoil.xTE,airfoil.yTE);
% find phi for boundary at boundary_distance away
% y_p = fzero( @(y) airfoil.phi_x_from_y(phiTE0,y,0).^2 + y.^2 - boundary_distance^2, boundary_distance );
% x_p = airfoil.phi_x_from_y(phiTE0,y_p,0);
% phiTEp = airfoil.phi_from_xy(airfoil.xTE,airfoil.yTE);
% tx = linspace(0,500,10001);
% ty = linspace(-500,500,10001);
% plot(tx,airfoil.psi_y_from_x(500*airfoil.chord,tx,-0.1),'g')
% plot(airfoil.phi_x_from_y(2.22,ty,0),ty,'b')
% plot(airfoil.phi_x_from_y(500*airfoil.chord,ty,0),ty,'b')


tic;
[x2,y2] = hyperbolic_C_grid_local_v2( x_airfoil, y_airfoil,                 ...
                                      n_wake_pts, jmax, ...
                                      boundary_distance=boundary_distance, ...
                                      wake_multiplier=wake_multiplier, ...
                                      AR_LE = AR_LE, ...
                                      AR_TE = AR_TE, ...
                                      scjmax = scjmax, ...
                                      mu = mu, ...
                                      muim = muim, ...
                                      alpham = alpham, ...
                                      fy_wake = fy_wake );
time = toc;
% hold on
% plot(x2,y2,'k');
% plot(x2.',y2.','k');
% axis equal
plot_edge_length_ratio(x2,y2,1,'jet',1);
xlim([-1,2])
ylim([-1.5,1.5])

plot_edge_length_ratio(x2,y2,2,'jet',1);
xlim([-1,2])
ylim([-1.5,1.5])

% L = edge_length_ratio_2D(x,y,dir);

hold off
%%
clf;
hold on
for j = 1:jmax
    tmp = sqrt(diff(x2(1:end-1,j)).^2+diff(y2(1:end-1,j)).^2) ./ sqrt(diff(x2(2:end,j)).^2+diff(y2(2:end,j)).^2);
    tmp(tmp<1) = 1./tmp(tmp<1);
    plot(tmp)
end

clf;
hold on
for i = 1:imax
    tmp = sqrt(diff(x2(i,1:end-1)).^2+diff(y2(i,1:end-1)).^2) ./ sqrt(diff(x2(i,2:end)).^2+diff(y2(i,2:end)).^2);
    tmp(tmp<1) = 1./tmp(tmp<1);
    plot(tmp)
end
hold off;


function hfig = plot_edge_length_ratio(x,y,dir,cmap,linewidth)
hfig = figure();
set(gca,'Color','k');
colormap(cmap);
L = edge_length_ratio_2D(x,y,dir);
% x = padarray(x,[1,1],nan,'post');
% y = padarray(y,[1,1],nan,'post');
% L = padarray(L,[1,1],nan,'post');
x1 = x; x1(end,:) = nan;
x2 = x; x2(:,end) = nan;
y1 = y; y1(end,:) = nan;
y2 = y; x2(:,end) = nan;

hold on;
patch(x,y,L,'edgecolor','flat','LineWidth',linewidth);
patch(x.',y.',L.','edgecolor','flat','LineWidth',linewidth);
axis equal
colorbar;
hold off
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