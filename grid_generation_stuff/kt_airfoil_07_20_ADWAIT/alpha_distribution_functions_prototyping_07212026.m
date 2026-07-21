%% 2D distribution functions (for alpham & scjmax) 07/21/2026
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


t = transfinite_interpolant_2D(G1y=[-1,-0.5,-0.5,-1], G3y=@(t)smooth_ramps(t,[0,1/3,15/32,17/32,2/3,1],[1.0000,0.999,1.0,1.0,0.999,1.0000]));

t.plot_grid(11);

function 


function f = bivariate_ramps(imax,jmax,i_vals,i_breaks,j_vals,j_breaks)
f_i = @(i) smooth_ramps((i-1)/(imax-1),i_breaks,i_vals);
f_j = @(j) smooth_ramps((j-1)/(jmax-1),j_breaks,j_vals);
f = @(i,j) f_i(i).*f_j(j);
end