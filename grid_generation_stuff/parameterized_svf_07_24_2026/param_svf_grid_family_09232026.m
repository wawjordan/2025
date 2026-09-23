%% svf AR param grid family (09/23/2026)
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
dstr = char(datetime('now',"Format",'uuuu-MM-dd''_''HH-mm-ss'));

folder = 'C:\Users\Will\Downloads\svf_for_transfer_AR_1000';
prefix='hb';
out_folder = fullfile(folder,'\grids\');
jobfmt  = ['_',prefix,'%0.4dx%0.4d'];

r_factor    = 2;
levels      = 1:6;
n_levels    = length(levels);

n_theta     = 97;
n_r         = 97;
AR_target   = 1000;

delta_s_max = 1.2;
theta_0     = pi/2;
theta_1     = 0;
r_0         = 2;
r_1         = 3;
INFO = generate_spacing_info( n_theta, n_r, theta_0, theta_1, r_0, r_1, AR_target, delta_s_max );

r = r_factor^(levels(end)-1);

% fine grid
imax = (INFO.n_theta-1)*r + 1;
jmax = (INFO.n_r    -1)*r + 1;
GRID = generate_grid( imax, jmax, INFO.theta_fun, INFO.r_fun );

bc_id_list  = [201,-200,201,-200];

Nfine = [GRID.imax,GRID.jmax];
for j = 1:n_levels
    s = r_factor^(levels(j)-1);
    fprintf('writing level %d of %d\n',j,n_levels)
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [out_folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\',prefix];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
    Ni(1) = 1;
    Ni(2) = (GRID.imax - 1)/s + 1;
    Nj(1) = 1;
    Nj(2) = (GRID.jmax - 1)/s + 1;
    N = [GRID2.imax,GRID2.jmax];
    idx_list{1} = [ Ni(1), Ni(2),  Nj(1), Nj(1) ];
    idx_list{2} = [ Ni(1), Ni(2),  Nj(2), Nj(2) ];
    idx_list{3} = [ Ni(1), Ni(1),  Nj(1), Nj(2) ];
    idx_list{4} = [ Ni(2), Ni(2),  Nj(1), Nj(2) ];
    additional_inputs{1} = [];
    additional_inputs{2} = [];
    additional_inputs{3} = [];
    additional_inputs{4} = [];
    SENSEI_BC_write(N,idx_list,bc_id_list,additional_inputs,[filename,'.bc'])
end

function GRID = generate_grid( n_theta, n_r, theta_fun, r_fun )
GRID = struct();
GRID.dim  = 2;
GRID.imax = n_theta;
GRID.jmax = n_r;
r     = r_fun(     linspace(0,1,n_r    ) );
theta = theta_fun( linspace(0,1,n_theta) );
[R,THETA] = ndgrid(r,theta);
GRID.x = R.*cos(THETA);
GRID.y = R.*sin(THETA);
end

function INFO = generate_spacing_info( n_theta, n_r, theta_0, theta_1, r_0, r_1, AR_target, delta_s_max )
INFO = struct();
INFO.n_theta   = n_theta;
INFO.n_r       = n_r;
INFO.theta_0   = theta_0;
INFO.theta_1   = theta_1;
INFO.r_0       = r_0;
INFO.r_1       = r_1;
INFO.AR_in     = AR_target;
INFO.dsmax_in  = delta_s_max;
refine = false;
[ INFO.dtheta_0, INFO.AR, INFO.theta_fun, INFO.r_fun, INFO.dsmax_out ] = optimize_spacing( n_theta, n_r, theta_0, theta_1, r_0, r_1, delta_s_max, AR_target, refine );
end

function [ dtheta_0, AR, theta_fun, r_fun, delta_s_max_out ] = optimize_spacing( n_theta, n_r, theta_0, theta_1, r_0, r_1, delta_s_max, AR_target, refine )

theta_total = abs(theta_1-theta_0);
eta         = r_1/r_0;
objfun = @(dt) ( AR_target - AR_mid(dt,theta_total,n_theta,refine) ).^2;

dt0 = theta_total/(n_theta-1);

A = [];
b = [];
lb = dt0/10000;
ub = dt0*2;
Aeq = [];
beq = [];
nonlcon = @(dt) dslecon(dt,delta_s_max,n_r,n_theta,eta,theta_total,refine);

% options = optimoptions("fmincon",OptimalityTolerance=1e-12,StepTolerance=1e-12);

dtheta_0 = fmincon(objfun,dt0,A,b,Aeq,beq,lb,ub,nonlcon);

[AR,deta0,delta_s_max_out] = AR_mid(dtheta_0,theta_total,n_theta,refine);

d0_t = dtheta_0/theta_total;
tau_theta = vinokur_two_sided_spacing_fcn(n_theta,d0_t,d0_t,refine);
theta_fun = @(t) (theta_1-theta_0)*tau_theta(t) + theta_0;

d0_r = log(deta0)/log(eta);
tau_r = vinokur_two_sided_spacing_fcn(n_r,d0_r,d0_r,refine);
r_fun = @(t) (r_1-r_0)*( eta.^tau_r(t) - 1 )/( eta - 1 )  + r_0;

function [ineqnonlin,eqnonlin] = dslecon(dtheta_0,dsmax,n_r,n_theta,eta,theta_total,refine)
    dsmax_ = max_stretching_rate(dtheta_0,n_r,n_theta,eta,theta_total,refine);
    ineqnonlin = dsmax_ - dsmax;
    eqnonlin = [];
end

function dsmax_out = max_stretching_rate(dtheta_0,n_r,n_theta,eta,theta_total,refine)
    [~,deta,dsmax_theta] = AR_mid(dtheta_0,theta_total,n_theta,refine);
    dsmax_r = get_radial_spacing_stretching_rate(eta,deta,n_r,refine);
    dsmax_out = max(dsmax_theta,dsmax_r);
end

end

function dsmax = get_radial_spacing_stretching_rate(eta,deta,n_r,refine)
d0 = log(deta)/log(eta);
[s0,~] = vinokur_two_sided_set_both_spacing(n_r,d0,d0,refine);
[delta,branch] = vinokur_get_delta_two_sided(s0,s0,refine);
[tau,~,~] = vinokur_two_sided_f(linspace(0,1,n_r),s0,s0,delta,branch);
xi = ( eta.^tau - 1 )/( eta - 1 );
s_rate  = diff(xi(2:end))./diff(xi(1:end-1));
dsmax = max(max(s_rate),1/min(s_rate));
end

function [AR,deta,dsmax] = AR_mid(dtheta_0,theta_total,n_theta,refine)
[dtheta_mid,dsmax] = arrayfun(@(dt0)get_dtheta_mid(dt0,theta_total,n_theta,refine),dtheta_0);
deta = (1+tan(dtheta_0/2))./(1-tan(dtheta_0/2));
AR   = tan(dtheta_mid/2) .* (deta+1)./(deta-1);
end

function [dtheta_mid,dsmax] = get_dtheta_mid(dtheta_0,theta_total,n_theta,refine)
d0 = dtheta_0/theta_total;
t0 = 0.5 - 1/(n_theta-1);
t1 = 0.5;
[s0,~] = vinokur_two_sided_set_both_spacing(n_theta,d0,d0,refine);
[delta,branch] = vinokur_get_delta_two_sided(s0,s0,refine);
[tau,~,~] = vinokur_two_sided_f([t0,t1],s0,s0,delta,branch);
if (nargout>1)
    t = vinokur_two_sided_f(linspace(0,1,n_theta),s0,s0,delta,branch);
    s_rate = diff(t(2:end))./diff(t(1:end-1));
    dsmax = max(max(s_rate),1/min(s_rate));
end
dtau = tau(2) - tau(1);
dtheta_mid = theta_total*dtau;
end

function P2D_grid_out(GRID,filename)
imax = GRID.imax;
jmax = GRID.jmax;
fid = fopen(filename,'w');
intfmt = '        %4d';
fltfmt = ' %-# 23.16E';
% fprintf(fid,[intfmt,'\n'],1); % this code is just for single block grids
fprintf(fid,[intfmt,intfmt,'\n'],imax,jmax);

count = 0;
for j = 1:jmax
    for i = 1:imax
        fprintf(fid,fltfmt,GRID.x(i,j));
        count = count + 1;
        if count == 2
            fprintf(fid,'\n');
            count = 0;
        end
    end
end

for j = 1:jmax
    for i = 1:imax
        fprintf(fid,fltfmt,GRID.y(i,j));
        count = count + 1;
        if count == 2
            fprintf(fid,'\n');
            count = 0;
        end
    end
end

fclose(fid);

end

function GRID = grid_subset_2D(GRID,small_ind)
GRID.x = GRID.x(small_ind{:});
GRID.y = GRID.y(small_ind{:});
GRID.imax = length(small_ind{1});
GRID.jmax = length(small_ind{2});
end