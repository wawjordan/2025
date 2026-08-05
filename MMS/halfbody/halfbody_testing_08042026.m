%% halfbody grid testing 08/04/2026
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

clc; clear; close all;

xstag = -10;
vinf   = 75;
length_scale = 1;
boundary_distance = 50;
mu = 0.4;
AR = 10;
% imax = 513;
% jmax = 129;
imax = 257;
jmax = 129;

n_ref = 0;


f = halfbody(xstag,vinf,length_scale);
GRID = f.extruded_grid(imax,jmax,boundary_distance,mu,AR);
hold on
plot(GRID.x(1:2^n_ref:end,1:2^n_ref:end),  GRID.y(1:2^n_ref:end,1:2^n_ref:end),  'k');
plot(GRID.x(1:2^n_ref:end,1:2^n_ref:end).',GRID.y(1:2^n_ref:end,1:2^n_ref:end).','k');
axis equal