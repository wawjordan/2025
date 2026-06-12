%% test script (06/04/2026)
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



xd =  [  0, 20,  40,  60,  80, 100 ];
N  =  [   70,   35,   50,   70,   110 ];

N = [1,cumsum(N)];
% N  =  [ 1, 21,  35,  61,  81,  90 ];
d  =  [ 0.01, 0.1,   0.01,   0.1,   1,   1 ];

i = 4;
xd = xd(1:i);
N  =  N(1:i);
d  =  d(1:i);

si = get_slopes(N,xd,d);

[x_,dx_,t_,xb,dxb,tb] = combine_curves(N,xd,d);
hold on
% plot(t_,x_,'k')
% plot(tb,xb,'r.')

plot(t_,dx_,'k')
plot(tb,dxb,'r.')

% plot(t_,dx_,'k')
plot(t_(1:end-1),diff(x_),'r')

hold off


function si = get_slopes_local(N_i,x_i,d_i)
% given a list of abscissae (ti) and spacings (di) compute a C1 stretching function
n_sub = numel(N_i) - 1;
N_i = N_i(:); % node count at xi (Ni(1)=1)
x_i = x_i(:); % location in physical space at Ni
d_i = d_i(:); % spacing at control point
si = zeros(n_sub,2);
si2 = zeros(n_sub,2);
refine = true;

L = diff(x_i);
dx = L./( diff(N_i) + 1);

d_i(1)  = d_i(1)/(x_i(2)-x_i(1));
for i = 2:n_sub+1
    d_i(i)  = d_i(i)/(x_i(i)-x_i(i-1));
end

% 1st subdomain
N = N_i(2) - N_i(1) + 1;
% L = ti(2) - ti(1);
d1 = d_i(1);
d2 = d_i(2);
[si(1,1),si(1,2)] = vinokur_two_sided_set_both_spacing(N,d1,d2,refine);
for i = 2:n_sub
% [s0,s1] = vinokur_two_sided_set_one_spacing(N,d0,s1,refine)
    N = N_i(i+1)-N_i(i) + 1;

    % rat = (N_i(i+1)-N_i(i))/(N_i(i)-N_i(i-1));
    d2 = d_i(i+1);
    si(i,1) = si(i-1,2);
    % s1 = si(i,1)/rat;
    s1 = si(i,1)/dx(i);
    [si(i,2),~] = vinokur_two_sided_set_one_spacing(N,d2,s1,refine);
end


%% now go in reverse
% last subdomain
N = N_i(n_sub+1) - N_i(n_sub) + 1;
d1 = d_i(n_sub+1);
d2 = d_i(n_sub);
[si2(n_sub,2),si2(n_sub,1)] = vinokur_two_sided_set_both_spacing(N,d1,d2,refine);
for i = n_sub-1:-1:1
% [s0,s1] = vinokur_two_sided_set_one_spacing(N,d0,s1,refine)
    N = N_i(i+1)-N_i(i) + 1;

    % rat = (N_i(i+1)-N_i(i))/(N_i(i)-N_i(i-1));
    d2 = d_i(i);
    si2(i,2) = si2(i+1,1);
    % s1 = si(i,1)/rat;
    s1 = si2(i,2)/dx(i);
    [si2(i,1),~] = vinokur_two_sided_set_one_spacing(N,d2,s1,refine);
end
% si = (si + si2)/2;

end



% function xi = stretch(t,Ni,ti,di)

% F1 = @(x) arrayfun( @(i) Fx(i).F(x), discretize(bound(x,y(1),y(end)),y) );