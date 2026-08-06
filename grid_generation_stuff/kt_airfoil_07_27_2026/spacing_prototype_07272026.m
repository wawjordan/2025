%% Apply stretching function to a given point distribution (07/27/2026)
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

% goal:
% given a nonuniform (but still monotonic) point distribution, and a
% desired target spacing at one or two ends of that distribution, calculate
% the stretching function that will give you this
N = 33;
spline_order = 5;
dx_ratio=0.1;
dx0=0.1;
x  = my_geomspace(N,0.0,xmax=20,dx0=dx0);
t1 = linspace(0,1,N);
xs = spapi( spline_order, t1, x );
dt0 = fzero( @(t)fnval(xs,t)-x(1)-dx_ratio*dx0,t1(2)*dx_ratio);

ft2 = vinokur_two_sided_spacing_fcn(N,1/(N-1),dt0,false);

t2 = ft2(t1);

x2 = fnval(xs,t2);

hold on;
plot(x,'k.-')
plot(x2,'r.-')



