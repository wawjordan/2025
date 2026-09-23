%% radius-aspect ratio function attempt 5 (09/22/2026)
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

theta_0   = pi/2;
theta_1   = 0;
r_0       = 2;
r_1       = 3;
n_theta   = 97;
n_r       = 97;
AR_target = 10000;
% n_theta   = 45;
% n_r       = 45;
% AR_target = 10;
dsmax     = 1.2;
refine    = false;

[dtheta_0,AR,f_theta,f_r,dsmax_out] = optimize_spacing(n_r,n_theta,r_0,r_1,theta_0,theta_1,dsmax,AR_target,refine);
% create grid
theta = f_theta( linspace(0,1,n_theta) );
r     = f_r(     linspace(0,1,n_r) );
[R,THETA] = ndgrid(r,theta);
X = R.*cos(THETA);
Y = R.*sin(THETA);

% calculate approximate aspect ratios
N = 8;
theta = f_theta( linspace(0,1,N*(n_theta-1)+1) );
r     = f_r(     linspace(0,1,N*(n_r-1)+1) );
eta_out    = r(2:end)./r(1:end-1);
dtheta_out = abs(diff(theta)).';
AR_out  = tan(dtheta_out/2).*( (eta_out+1)./(eta_out-1) );


% plot
tiledlayout(1,2)
nexttile
hold on
plot(X,Y,'r')
plot(X.',Y.','r')
axis equal

nexttile
AR_plot = log10(max(AR_out,1./AR_out));
AR_plot = padarray(AR_plot,[1,1],'post');
pcolor(AR_plot,"FaceColor","flat","EdgeColor","none");
colorbar


function [dtheta_0,AR,f_theta,f_r,dsmax_out] = optimize_spacing(n_r,n_theta,r0,r1,theta0,theta1,dsmax,AR_target,refine)

theta_total = abs(theta1-theta0);
eta         = r1/r0;
objfun = @(dt) ( AR_target - AR_mid(dt,theta_total,n_theta,refine) ).^2;

dt0 = theta_total/(n_theta-1);

A = [];
b = [];
lb = dt0/10000;
ub = dt0*2;
Aeq = [];
beq = [];
nonlcon = @(dt) dslecon(dt,dsmax,n_r,n_theta,eta,theta_total,refine);

% options = optimoptions("fmincon",OptimalityTolerance=1e-12,StepTolerance=1e-12);

dtheta_0 = fmincon(objfun,dt0,A,b,Aeq,beq,lb,ub,nonlcon);

[AR,deta0,dsmax_out] = AR_mid(dtheta_0,theta_total,n_theta,refine);

d0_t = dtheta_0/theta_total;
tau_theta = vinokur_two_sided_spacing_fcn(n_theta,d0_t,d0_t,refine);
f_theta = @(t) (theta1-theta0)*tau_theta(t) + theta0;

d0_r = log(deta0)/log(eta);
tau_r = vinokur_two_sided_spacing_fcn(n_r,d0_r,d0_r,refine);
f_r = @(t) (r1-r0)*( eta.^tau_r(t) - 1 )/( eta - 1 )  + r0;


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