%% testing for my_geomspace_w_refine (02/18/2026)
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


xmin = 0;
xmax_in = 100;
dx0_in = 0.00001;
dx1_in = 10;
r_in   = 1.00001;
N = 1;

hold on
% case 1
[x_old,dx0_old,dx1_old] = my_geomspace(N,xmin,xmax=xmax_in,r=r_in);
[x_new,dx0_new,dx1_new,~] = my_geomspace_w_refine(N,1,xmin,xmax=xmax_in,r=r_in);
plot(abs(x_old-x_new));
disp(dx0_old-dx0_new)
disp(dx1_old-dx1_new)


% case 2
[x_old,xmax_old,dx1_old] = my_geomspace(N,xmin,dx0=dx0_in,r=r_in);
[x_new,xmax_new,dx1_new,~] = my_geomspace_w_refine(N,1,xmin,dx0=dx0_in,r=r_in);
plot(abs(x_old-x_new));
disp(xmax_old-xmax_new)
disp(dx1_old-dx1_new)

% case 3
[x_old,xmax_old,r_old] = my_geomspace(N,xmin,dx0=dx0_in,dx1=dx1_in);
[x_new,xmax_new,r_new,~] = my_geomspace_w_refine(N,1,xmin,dx0=dx0_in,dx1=dx1_in);
plot(abs(x_old-x_new));
disp(xmax_old-xmax_new)
disp(r_old-r_new)

% case 4
[x_old,r_old,dx1_old] = my_geomspace(N,xmin,xmax=xmax_in,dx0=dx0_in);
[x_new,r_new,dx1_new,~] = my_geomspace_w_refine(N,1,xmin,xmax=xmax_in,dx0=dx0_in);
plot(abs(x_old-x_new));
disp(r_old-r_new)
disp(dx1_old-dx1_new)