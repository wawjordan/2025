%% prototype script for slip-wall BC for variational rec. (06/22/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 2;
neq = 4;
syms n [d 1]
syms t [d 1]

M1 = [ diag(n(1)*ones(neq,1)),diag(n(2)*ones(neq,1)); ...
       diag(t(1)*ones(neq,1)),diag(t(2)*ones(neq,1)) ];

M2 = diag([((-1).^(1:neq).');((-1).^(0:neq-1).')]) * M1;

A = M2\M1;
pretty(A)
% subs(A,{n(1),n(2),t(1),t(2)},{1,0,0,1})

