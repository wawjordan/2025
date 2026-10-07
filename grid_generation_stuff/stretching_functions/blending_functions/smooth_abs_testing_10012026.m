%% messing with smooth abs (10/01/2026)
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


t = linspace(0,1,1001);
a = 1/100;
N = 33;
[d,a_out] = vinokur_set_mid_spacing(N,a);

f = vinokur_two_sided_spacing_fcn( N, d, d, false);
plot(t,f(t))

[f0,df0,ddf0] = my_abs_fcn();

k = 1;
[f1,df1,ddf1] = smooth_abs_fcn(k);

k = 10;
[f2,df2,ddf2] = smooth_abs_fcn(k);

tiledlayout(1,3)

nexttile
hold on
fplot(@(x)f0(x),[-1,1],'k')
fplot(@(x)f1(x),[-1,1],'r')
fplot(@(x)f2(x),[-1,1],'b')

nexttile
hold on
fplot(@(x)df0(x),[-1,1],'k')
fplot(@(x)df1(x),[-1,1],'r')
fplot(@(x)df2(x),[-1,1],'b')

nexttile
hold on
fplot(@(x)ddf0(x),[-1,1],'k')
fplot(@(x)ddf1(x),[-1,1],'r')
fplot(@(x)ddf2(x),[-1,1],'b')


function [d,a_out] = vinokur_set_mid_spacing(N,a)
dt = 1/(N-1);
if (a>0.5)
    error('a too big')
end
% initial guess for slope
s0 = dt;
options = optimset();
s = fzero(@(s)obj_fun(s,dt,a),s0,options);
a_out = obj_fun(s,dt,a) + a;

[delta_,branch_] = vinokur_get_delta_two_sided(s,s,false);
[xi_,~] = vinokur_two_sided_f([0,dt],s,s,delta_,branch_);
d = xi_(2)-xi_(1);

function e = obj_fun(s,dt,a_target)
[delta,branch] = vinokur_get_delta_two_sided(s,s,false);
[xi,~] = vinokur_two_sided_f([0.5-dt,0.5],s,s,delta,branch);
a_tmp = xi(2)-xi(1);
e = a_tmp - a_target;
end
end