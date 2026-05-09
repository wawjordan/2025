%% Ringleb Flow prototyping - from "I do like CFD"
% modifications: 04/06/2026
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

N_theta = 61;
N_psi   = 31;

kmax = 0.8;
kmin = 0.3;
qbar = 0.2;

psi_min = 1/kmax;
psi_max = 1/kmin;
V_min   = qbar;


% subsonic
psi_min = 2;
psi_max = 6;
V_min   = 0.12;

% supersonic
% psi_min = 0.6;
% psi_max = 0.8;
% V_min   = 1.2;

% [X1,Y1,~,~] = basic_grid_1_full(N_theta,N_psi,psi_min,psi_max,V_min);

psi_sp1 = @(xmin,xmax,N) my_geomspace( N, xmin, xmax=xmax, r=1.08);
phi_sp1 = @(xmin,xmax,N) linspace(xmin,xmax,N);
phi_sp2 = @(xmin,xmax,N) sin_space(xmin,xmax,N);
[X1,Y1,~,~] = basic_grid_2_full(N_theta,N_psi,psi_min,psi_max,V_min,psi_sp1,phi_sp1);
[X2,Y2,~,~] = basic_grid_2_full(N_theta,N_psi,psi_min,psi_max,V_min,psi_sp1,phi_sp2);

hold on; axis equal;
plot(-X1,Y1,'r'); plot(-X1.',Y1.','r');
plot(-X2,Y2,'k'); plot(-X2.',Y2.','k');


function [X,Y,psi,V] = basic_grid_1_full(N_theta,N_psi,psi_min,psi_max,V_min)

% psi_0 = linspace( psi_min, psi_max, N_psi );
psi_0 = my_geomspace( N_psi, psi_min, xmax=psi_max, r=1.08);
V   = zeros(N_theta,N_psi);
psi = zeros(N_theta,N_psi);

for i = 1:N_theta
    psi(i,:) = psi_0;
end

for j = 1:N_psi
    theta_min = asin(psi_0(j)*V_min);
    theta_max = pi - asin(psi_0(j)*V_min);
    theta = linspace(theta_min,theta_max,N_theta);
    V(:,j) = sin(theta)/psi_0(j);
end

if mod(N_theta,2) == 1
    V((N_theta+1)/2,:) = 1.0./psi_0;
    V((N_theta+1)/2+1:N_theta,:) = -V((N_theta+1)/2+1:N_theta,:);
else
    V(N_theta/2+1:N_theta,:) = -V(N_theta/2+1:N_theta,:);
end

X = x_fun(psi,V);
Y = y_fun(psi,V);

end

function x = piecewise(xmin,xmax,n)
x = linspace(0,1,n);

x = exp(1)*(x(1:end-1)-1).*exp(1./(x(1:end-1)-1))+1;
x(end) = 1;
x = 2*x -1;
x = xmin + (xmax-xmin)*0.5*x;
end

function x = sin_space(xmin,xmax,n)
% x = xmin + (xmax-xmin)*0.5*( sin(linspace(-1,1,n)*pi/2) + 1);

% x = xmin + (xmax-xmin)*0.5*( ( sin(linspace(-1,0,n)*pi/4) / sin(pi/4) ).^3 + 1);
tmp = linspace(0,1,n);
x = xmin + (xmax-xmin)*0.5*( sin(tmp*pi/2)+ sin(tmp*pi/2).^7 );
end

function [X,Y,psi,V] = basic_grid_2_full(N_theta,N_psi,psi_min,psi_max,V_min, psi_space, phi_space )
psi_0 = psi_space( psi_min, psi_max, N_psi );
V   = zeros(N_theta,N_psi);
psi = zeros(N_theta,N_psi);

for i = 1:N_theta
    psi(i,:) = psi_0;
end

for j = 1:N_psi
    theta_min = asin(V_min*psi_0(j));
    phi_start = 1./(psi_0(j)*V_min)  *cos(theta_min);
    phi_end   = -phi_start;
    phi1 = linspace(phi_start,phi_end,N_theta);
    phi2 = piecewise(phi_start,phi_end,N_theta);
    phi = phi_space(phi_start,phi_end,N_theta);
    for i = 1:ceil((N_theta+1)/2)
        k = 1/psi_0(j);
        V(i,j) = sqrt(k^2/(1+phi(i)^2));
    end
end

% if mod(N_theta,2) == 1
%     V((N_theta+1)/2,:) = 1.0./psi_0;
%     V((N_theta+1)/2+1:N_theta,:) = -V((N_theta+1)/2+1:N_theta,:);
% else
%     V(N_theta/2+1:N_theta,:) = -V(N_theta/2+1:N_theta,:);
% end

X = x_fun(psi,V);
Y = y_fun(psi,V);

end


function b = b_value(V)
b = sqrt(1.0 - 0.2*V.^2);
end

function rho = density(V)
b = b_value(V);
rho = b.^5;
end

function psi = psi_fun(xp,V)
    b   = b_value(V);
    rho = b.^5;
    L   = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) ...
        - 0.5*log((1.0+b)./(1.0-b));
    psi = sqrt(0.5./(V.^2) - rho.*(xp-0.5*L));
end

function L = L_value(V)
b = b_value(V);
L = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) - 0.5*log((1.0+b)./(1.0-b)); 
end

function x = x_fun(psi,V)
rho = density(V);
L   = L_value(V);
x = 1.0./rho .* (0.5./(V.^2) - psi.^2) + 0.5*L; 
end

function y = y_fun(psi,V)
rho = density(V);
rad = 1.0 - (V.^2).*(psi.^2);
rad = (rad>=0).*rad + (rad<0)*0;
y = psi./(rho.*V) .* sqrt(rad);
end


function state = exact_soln_nondim(xp,V)
    b   = b_value(V);
    rho = b.^5;
    p   = b.^7;
    L   = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) ...
        - 0.5*log((1.0+b)./(1.0-b));
    psi   = sqrt(0.5./V.^2 - rho.*(xp-0.5*L));
    arg = psi*V;
    arg = arg*(abs(arg)<=1) + 1*(arg>1) - 1*(arg<-1);
    theta = asin(arg);
    u     = -V*cos(theta);
    v     = -V*sin(theta);
    state = 0.0;
    state(1) = rho;
    state(2) = u;
    state(3) = v;
    state(5) = p;
end


function S = get_soln(X,V)
[M,N] = size(X);
S = zeros(M,N,5);
for j = 1:N
    for i = 1:M
        S(i,j,:) = exact_soln_nondim(X(i,j),V(i,j));
    end
end
end




function V = find_solution_V(xp,yp,V0,maxiter,tol)
% initial guess
% Vm1 = 0.5*(V_min + 1/psi_min);
V = Vguess(xp,yp,V0);
err0 = abs(V-V0);
V0 = V;
for i = 1:maxiter
    V = Vguess(xp,yp,V0);
    R = abs(V-V0)/err0;
    if R < tol
        break
    end
    V0 = V;
end

if (R > tol)
    warning('tol not met in find_solution_V')
end

end

function V = Vguess(xp,yp,V)
    rho = density(V);
    L   = L_value(V);
    den = 2.0*rho*sqrt((xp-0.5*L)^2 + yp^2);
    V = sqrt(1.0/den);
end