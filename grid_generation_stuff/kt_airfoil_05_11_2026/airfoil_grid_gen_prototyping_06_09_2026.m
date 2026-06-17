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

N        = 257;
Ntrans   = 17;

AR_LE    = 1;
AR_TE    = 1;
% delta_LE = 0.0005*AR_LE;
% delta_TE = 0.0005*AR_TE;
delta_LE = 0.001*AR_LE;
delta_TE = 0.001*AR_TE;

[x_airfoil,y_airfoil,F1,dF1,ddF1] = get_airfoil_coordinates(airfoil,N,Ntrans,delta_LE,delta_TE);

hold on
plot(x_airfoil,y_airfoil,'r.-');
axis equal
hold off

boundary_distance = 500;
jmax              = 129;             % Number of Points Off Body
n_wake_pts        = 129;             % Number of Points in Wake
scjmax            = 0.9999;           % Scaling Factor
mu                = 0.1;             % 4th order explicit smoothing factor
muim              = 0.5;             % Implicit smoothing factor
alpham            = 0.5;            % Alpha scheme integration factor
% alpham            = 0.5;            % Alpha scheme integration factor


[x2,y2] = hyperbolic_C_grid_local_v2( x_airfoil, y_airfoil,                 ...
                                 boundary_distance, jmax, ...
                                 AR_LE, AR_TE, scjmax, ...
                                 mu, muim, alpham, ...
                                 n_wake_pts );

hold on
plot(x2,y2,'k');
plot(x2.',y2.','k');
axis equal

hold off
%%
% clf;
% hold on
% for j = 1:jmax
%     tmp = sqrt(diff(x2(1:end-1,j)).^2+diff(y2(1:end-1,j)).^2) ./ sqrt(diff(x2(2:end,j)).^2+diff(y2(2:end,j)).^2);
%     tmp(tmp<1) = 1./tmp(tmp<1);
%     plot(tmp)
% end
% imax = size(x2,1);
% tmp = zeros(imax,jmax-1);
% for i = 1:imax
%     tmp(i,:) = sqrt(gradient(x2(i,1:end-1)).^2+gradient(y2(i,1:end-1)).^2) ./ sqrt(gradient(x2(i,2:end)).^2+gradient(y2(i,2:end)).^2);
% end
% tmp(tmp<1) = 1./tmp(tmp<1);
% contourf(tmp);

function [x,y,F1,dF1,ddF1] = get_airfoil_coordinates(airfoil,N,Ntrans,delta_LE,delta_TE)
% get parametric location of maximum curvature on the airfoil
% which will approximately be the leading edge
tLE = airfoil.thetaCmax/(2*pi);

% get normalized arc length to point of max curvature
t0 = airfoil.airfoil_arc_length(0,tLE)/airfoil.airfoil_arc_length(0,1);

% nondimensionalize spacings
cxL   = airfoil.chord/airfoil.airfoil_arc_length(0,1);
dLE = delta_LE*cxL;
dTE = delta_TE*cxL;

% get buffer distance over which to blend the two functions
offset = ( (Ntrans-1)/(N-1) )*cxL;

% generate asymmetric stretching function for airfoil surface
[F1,dF1,ddF1] = hermite_blend_2_vinokur_asym(N,t0,dLE,dTE,offset,true);
% % [f,df,ddf] = hermite_blend_2_vinokur_asym(N,t0,d0,d1,off,refine)
[x,y] = airfoil.output_airfoil_coords1(N,F1);

% flip to match the clockwise convention
x = flip(x);
y = flip(y);
end

function airfoil = local_airfoil_generator()
epsilon = 0.1;
kappa   = 0.0;
tau     = 0.0;
vinf    = 75.0;
rhoinf  = 1.0;
pinf    = 100000.0;
gamma   = 1.4;
alpha   = 10; % (degrees)
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