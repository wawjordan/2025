%% radius-aspect ratio function attempt 4 (09/21/2026)
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
n_theta = 129;
n_r     = 145;
x0 = 2;
x1 = 3;

eta = x1/x0;
dtheta = theta/(n_theta-1);

AR0 = 1000;
AR1 = 1000;

eta0 = (AR0 + tan(dtheta/2))/(AR0 - tan(dtheta/2));
eta1 = (AR1 + tan(dtheta/2))/(AR1 - tan(dtheta/2));

delta0 = log(eta0)/log(eta);
delta1 = log(eta1)/log(eta);


% delta1 = 1 - log(eta/eta1)/log(eta);
% tau = vinokur_two_sided_spacing_fcn( n_r, delta0, delta1,true);

% flipped !!!??
tau = vinokur_two_sided_spacing_fcn( n_r, delta1, delta0,true);


% exponent is defined by the base number of radial nodes
beta = n_r-1;

% solve for the coefficient
alpha = eta^(1/beta);

% create function for normalized spacing
xi = @(t) ( alpha.^(beta*tau(t)) - 1 )/( eta - 1 );

% dimensioned spacing
x = @(t) (x1-x0)*xi(t) + x0;


% create grid
t = linspace(0,1,n_r);
r = x(t);
th = linspace(theta,0,n_theta);
[R,THETA] = ndgrid(r,th);
X = R.*cos(THETA);
Y = R.*sin(THETA);

% calculate approximate aspect ratios
eta_out = r(2:end)./r(1:end-1);
AR_out1 = ( (eta_out-1)./(eta_out+1) ) .* (2./dtheta);
AR_out2 = tan(dtheta/2)*( (eta_out+1)./(eta_out-1) );


% plot
tiledlayout(1,2)
nexttile
hold on
plot(X,Y,'r')
plot(X.',Y.','r')
axis equal

nexttile
hold on
plot(AR_out1)
plot(AR_out2)