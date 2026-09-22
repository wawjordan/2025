%% radius-aspect ratio function attempt 3 (09/20/2026)
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

theta = pi/2;
% n_theta = 101;
% n_r     = 17;
% n_theta = 95;
% n_r     = 25;
n_theta = 129;
n_r     = 33;
% n_theta = 2049;
% n_r     = 513;
x0 = 2;
x1 = 3;

eta = x1/x0;

dtheta = theta/(n_theta-1);

dr0    = x0*dtheta;
dr1    = x1*dtheta;

% delta  = 0.5*(dr0+dr1)/(x1-x0);


% set these to evaluate at more points on the existing spacing functions
% AR_r = 1/2;
% AR_t = 1/2;
AR_r = 1;
AR_t = 1;

% exponent is defined by the base number of radial nodes
beta = n_r-1;

% solve for the coefficient
alpha = eta^(1/beta);

% estimate normalized spacing for first d0 and d1 in the stretching
% function, and average
delta0 = log(dr0*(eta-1)+1)/(beta*log(alpha));
delta1 = 1-log((1-dr1)*(eta-1)+1)/(beta*log(alpha));
delta  = 0.5*(delta0+delta1);

% create (symmetric) stretching function
tau = vinokur_two_sided_spacing_fcn( n_r, delta, delta,false);

% create function for normalized spacing
xi = @(t) ( alpha.^(beta*tau(t)) - 1 )/( eta - 1 );

% dimensioned spacing
x = @(t) (x1-x0)*xi(t) + x0;

% compute new node counts
n_r2     = AR_r*(n_r-1)+1;
n_theta2 = AR_t*(n_theta-1)+1;

% create grid
t = linspace(0,1,n_r2);
r = x(t);
th = linspace(theta,0,n_theta2);
[R,THETA] = ndgrid(r,th);
X = R.*cos(THETA);
Y = R.*sin(THETA);

% calculate approximate aspect ratios
dtheta2 = theta/(n_theta2-1);
alpha_out = r(2:end)./r(1:end-1);
AR_out = ( (alpha_out-1)./(alpha_out+1) ) .* (2./dtheta2);

% plot
tiledlayout(1,2)
nexttile
hold on
plot(X,Y,'r')
plot(X.',Y.','r')
axis equal

nexttile
plot(AR_out)