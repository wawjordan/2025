%% testing reparam_complex_curve (07/16/2026)
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
in.epsilon = 0.1;
in.kappa   = 0.0;
in.tau     = 0.1;
airfoil = kt_airfoil( in.epsilon, in.kappa, in.tau );
N = 33;
n_iter = 100;
tmin = 0;
tmax = 2*pi;
mu = 0.5;

% fs = @(t) t;
fs = vinokur_two_sided_spacing_fcn(N,0.1/(N-1),0.1/(N-1),false);
fz = @(theta) airfoil.z_from_theta(theta);
t0 = linspace(0,1,N);
t_exact = airfoil.arc_length_param2(t0,fs);
L_exact = airfoil.airfoil_arc_length(0,1);

% hold on
% plot(t0,t_exact,'k');
% plot(t0,t_approx,'r');

err_t = zeros(n_iter,3);
err_L = zeros(n_iter,1);
for i = 1:n_iter
    [t_approx,L_approx] = reparam_complex_curve_local(fz,fs,N,tmin,tmax,i,mu);
    t_approx = t_approx(:);
    err = abs(t_approx(2:end-1)-t_exact(2:end-1))/L_exact;
    err_t(i,1) = sum(err)/(N-2);
    err_t(i,2) = sqrt( sum(err.^2)/(N-2) );
    err_t(i,3) = max(err);
    err_L(i,1) = abs(L_approx-L_exact)/L_exact;
end

hold on
plot(err_L(:,1),'k.-')
plot(err_t(:,1),'r.-')
plot(err_t(:,2),'b.-')
plot(err_t(:,3),'g.-')
set(gca,'YScale','log')
legend({'length error','L_1 t error','L_2 t error','L_inf t error'})

for i = 1:n_iter
    [t_approx,L_approx] = reparam_complex_curve_local(fz,fs,N,tmin,tmax,i,1);
    t_approx = t_approx(:);
    err = abs(t_approx(2:end-1)-t_exact(2:end-1))/L_exact;
    err_t(i,1) = sum(err)/(N-2);
    err_t(i,2) = sqrt( sum(err.^2)/(N-2) );
    err_t(i,3) = max(err);
    err_L(i,1) = abs(L_approx-L_exact)/L_exact;
end

hold on
plot(err_L(:,1),'ko--')
plot(err_t(:,1),'ro--')
plot(err_t(:,2),'bo--')
plot(err_t(:,3),'go--')
set(gca,'YScale','log')


function [t,L] = reparam_complex_curve_local(f1,f2,N,tmin,tmax,n_iter,mu)
% creates a set of points parameterized in (approximate) arc length space
% of f1(tmin,tmax) with spacing from f2(0,1)
t0 = linspace(0,1,N);
t = t0;
for i = 1:n_iter
    z = f1( tmin + (tmax-tmin)*t ); z = z(:).';
    points = [real(z);imag(z)];
    % tc = [ 0; cumsum( sqrt( sum( abs(points(:,2:end) - points(:,1:end-1)).^2, 1) ) ).'];
    tc = [ 0; cumsum( sqrt( sum( abs(points(:,2:end) - points(:,1:end-1)).^2, 1).^mu ) ).'];
    L = tc(end);
    t = interp1(tc/L,t,f2(t0));
end
end