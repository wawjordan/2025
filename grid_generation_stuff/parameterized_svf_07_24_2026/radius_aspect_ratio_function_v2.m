%% radius-aspect ratio function attempt 2 (09/18/2026)
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
x0 = 2;
x1 = 3;

eta = x1/x0;

dtheta = theta/(n_theta-1);

% AR = @(t) 1-0.2*sin(pi*t);
AR = 1;
y = (dtheta/2)*AR;
alpha = (1 + y)./(1 - y);
beta = log( eta )/log( alpha );

% tau = @(t) log(1+t*(exp(1)-1));
% tau = @(t) ( log( (eta-1)*t + 1 )/(beta*log(alpha)) );

tau = @(t) (t - cos(pi*t)+1)/3;

tau = vinokur_two_sided_spacing_fcn( n_theta, 0.001, 0.001,false);

xi = @(t) ( alpha.^(beta*tau(t)) - 1 )/( eta - 1 );

x = @(t) (x1-x0)*xi(t) + x0;

% fplot(@(t)x(t),[0,1])

n_r = floor(beta)+1;

t = linspace(0,1,n_r);

r = x(t);

alpha = r(2:end)./r(1:end-1);
AR_out = ( (alpha-1)./(alpha+1) ) .* (2./dtheta);
% hold on
% plot(AR_out)
% plot(AR_out./AR(t(1:end-1)))

[R,THETA] = ndgrid(r,linspace(0,theta,n_theta));

X = R.*cos(THETA);
Y = R.*sin(THETA);

hold on
plot(X,Y,'r')
plot(X.',Y.','r')
axis equal