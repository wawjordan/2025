%% Smoothmax Testing (07/01/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
grid_file = fullfile(parent_dir,'kt.grd');
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

f1   = @(x) x.^2;
df1  = @(x) 2*x;
ddf1 = @(x) 2*ones(size(x));

f2   = @(x)  sin(x);
df2  = @(x)  cos(x);
ddf2 = @(x) -sin(x);

e = 0.001;
[f,df,ddf] = smoothmax_smu_fcn(f1,f2,df1,df2,ddf1,ddf2,e);

hold on
fplot(f1,[-1,1],'r')
fplot(f2,[-1,1],'b')
fplot(f,[-1,1],'g')
hold off
% [f,df,ddf] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);


function [f,df,ddf] = smoothmax_smu_fcn(f1,f2,df1,df2,ddf1,ddf2,e)
f   = @(x)   sm_f(x,f1,f2,df1,df2,ddf1,ddf2,e);
df  = @(x)  dsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e);
ddf = @(x) ddsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e);
end

function f = sm_f(x,f1,f2,df1,df2,ddf1,ddf2,e)
[f,~,~] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);
end

function df = dsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e)
[~,df,~] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);
end

function ddf = ddsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e)
[~,~,ddf] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);
end