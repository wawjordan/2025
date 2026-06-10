%% hermite_blended_piecewise_linear testing (06/10/2026)
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
% x = [0,0.2,1.0,1.8,2.6];
% y = [0,0.5,0.7,0.8,2.0];
% dm = [0.10,0.10,0.10];
% dp = [0.10,0.10,0.10];
% f = hermite_blended_piecewise_linear(x,y,dm,dp);
% fplot(f,[-1,3])


clf;
o   = 10;
N1 = 257;
N   = [(N1-1)/2+1-o,(N1-1)/2+1+o];
xi1 = [0,cumsum(N-1)/sum(N(:)-1)];
xi2 = [0,0.3,1.0];
d   = diff(xi2)./(N-1);

f = hermite_blended_piecewise_linear(xi1,xi2,10*d(1),10*d(2));
fplot(f,[0,1])