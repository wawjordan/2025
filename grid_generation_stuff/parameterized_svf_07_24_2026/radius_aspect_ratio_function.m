%% radius-aspect ratio function (09/17/2026)
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
AR = @(t) 1+zeros(size(t));
AR = @(t) 1 - 0.2*t;

y = @(t) (dtheta/2)*AR(t);

alpha = @(t) (1 + y(t))./(1 - y(t));

beta = log( eta )/log( alpha(1) );

xi = @(t) ( alpha(t).^(beta*t) - 1 )/( eta - 1 );

x = @(t) (x1-x0)*xi(t) + x0;

% fplot(@(t)x(t),[0,1])

n_r = floor(beta)+1;

t = linspace(0,1,n_r);

r = x(t);

[R,THETA] = ndgrid(r,linspace(0,theta,n_theta));

X = R.*cos(THETA);
Y = R.*sin(THETA);

hold on
plot(X,Y,'r')
plot(X.',Y.','r')
axis equal