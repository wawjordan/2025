%% test script: blending functions from Ashby's Thesis (06/09/2026)
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


% N = 33;
% Ni = 10;
% A = 0.9;
% t = linspace(0,N-1,1001);
% 
% f1 = @(t) exp(-t.^2);
% f2 = @(t) 0.5 + 0.5*sin(t);
% f3 = @(t) sqrt(t/(N-1)).*exp(-t./(N-1));
% 
% b1 = @(t) rbf_interior_1D(t,1,33,Ni,A);
% b2 = @(t) 1 - rbf_interior_1D(t,1,33,Ni,A);
% b3 = @(t) smooth_transition(t,0,10,11,15);


t = linspace(0,1,1001);
% val = smooth_bump_scale(x,xloc,w,wm,wp,scale_factor)
% f = combine_functions(f,w)
b1 = @(t) smooth_bump_scale(t,0.25,0.0,0.05,0.05,0.5);
b2 = @(t) smooth_bump_scale(t,0.75,0.0,0.05,0.05,0.5);
f = @(t) smooth_transition(t,0,0.25,0.75,1);
% b = combine_functions({b1,b2},{@(t)1-b1(t),@(t)1-b2(t)});
b = combine_functions({@(t)1-b1(t),@(t)1-b2(t)},{b1,b2});

% f1 = combine_functions({f},{b});
% f1 = combine_functions({f,f},{b1,b2});
f1 = @(t) f(t) + b1(t) + b2(t)-2;
hold on
plot(t,b1(t),'r--')
plot(t,b2(t),'b--')
plot(t,b(t),'g')
plot(t,f(t),'k')
plot(t,f1(t),'m')


function val = rbf_interior_1D(xi,xi_max,Nmax,Ni,A)
rmax = 0.5*(xi_max+1);
% rd_  = sqrt( ((Ni-1)/(Nmax-1)-0.5*(xi_max+1)).^2 );
rd_  = sqrt( (Ni-0.5*(Nmax+1)).^2 );
d    = -log(A)*(rmax/rd_)^2;
rd   = sqrt( (xi-0.5*(xi_max+1)).^2 );
val  = exp(-d*(rd/rmax).^2);
end