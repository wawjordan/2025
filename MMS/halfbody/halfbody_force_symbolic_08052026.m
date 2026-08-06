%% force on the halfbody 08/05/2026
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

syms t real % theta
syms v real % vinf
syms s real % sigma
syms p real % pinf
syms d real % rhoinf
sigma = 1;
vinf  = 75;
pinf  = 100000;
rhoinf = 1;

s_ = s/(2*pi);

r  = (s_/v)*(pi-t)/sin(t);

drdt = diff(r,t);
% simplify(drdt==-s_*((pi-t)*cos(t) + sin(t))/(sin(t)^2));

uvel  = v + s_*cos(t)/r;
vvel  =   - s_*sin(t)/r;
v2    =  uvel^2 + vvel^2;
nr    =  r;
nt    = -drdt;
nmag  = sqrt(r^2+drdt^2);
exr   =  cos(t);
ext   = -sin(t);

% exnhat = (nr*exr + nt*ext)/nmag;
% nmag gets canceled out with the differential line element
exnhat = (nr*exr + nt*ext);

integrand_ = (p-p-(1/2)*d*(v^2 - v2))*simplify(exnhat);
integrand_ = simplify(subs(integrand_,{s,v,d,p},{sigma,vinf,rhoinf,pinf}),100);

intfun = matlabFunction(integrand_);

val = 2*vpaintegral(integrand_,0,pi);
val2 = 2*integral(intfun,0,pi);

% p = pinf + 0.5*rhoinf*()
% Cx


fplot(intfun,[0,2*pi])
