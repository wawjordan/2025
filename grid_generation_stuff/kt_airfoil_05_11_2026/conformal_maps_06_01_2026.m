%% conformal mapping shenanigans (06/01/2026)
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

inputs = struct();
inputs.epsilon = 0.1;
inputs.kappa   = 0.0;
% inputs.kappa   = 0.1;
inputs.tau     = 0.0;
inputs.vinf    = 75.0;
inputs.rhoinf  = 1.0;
inputs.pinf    = 100000.0;
inputs.gamma   = 1.4;

inputs.alpha   = 5; % (degrees)
inputs.nskip   = 4;
inputs.rho_ref = 1.0;
inputs.p_ref   = 100000.0;
inputs.a_ref   = sqrt(inputs.gamma*inputs.p_ref/inputs.rho_ref);

% nondimensionalize inputs
inputs.vinf   = inputs.vinf/inputs.a_ref;
inputs.rhoinf = inputs.rhoinf/inputs.rho_ref;
inputs.pinf   = inputs.pinf/(inputs.rho_ref*inputs.a_ref^2);

airfoil        = kt_airfoil( inputs.epsilon, inputs.kappa, inputs.tau );
airfoil.vinf   = inputs.vinf;
airfoil.rhoinf = inputs.rhoinf;
airfoil.pinf   = inputs.pinf;
airfoil        = airfoil.set_alpha(inputs.alpha);

Ni = 50;
Nj = 50;

%% Define the initial grid
x0 = linspace(-pi,pi,Ni);
y0 = linspace(-pi,pi,Nj);
% x0 = linspace(-pi/2,pi/2,Ni);
% y0 = linspace(-pi/2,pi/2,Nj);
[x1,y1] = ndgrid(x0,y0);
z1 = x1 + 1i*y1;

%% Apply first mapping
% z2 = (z1-1)./(z1+1);
z2 = exp(z1);
z2 = 1./(4*z2) + 3*z2/4;
% z2 = exp(z1);
x2 = real(z2);
y2 = imag(z2);

% z3 = (1i*(1-z1)./(z1+1));
% z3 = log(z2);
% z3 = (1i-z1)./(z1+1i);
% z3 = airfoil.zeta_to_z(z2);
z3 = slit_disc(z2);
x3 = real(z3);
y3 = imag(z3);

figure
subplot(1,3,1)
hold on
plot(x1,y1,'r')
plot(x1.',y1.','b')
plot(x1(:,1).',y1(:,1).','g')
plot(x1(:,end).',y1(:,end).','m--')

subplot(1,3,2)
hold on
plot(x2,y2,'r')
plot(x2.',y2.','b')
plot(x2(:,1).',y2(:,1).','g')
plot(x2(:,end).',y2(:,end).','m--')


subplot(1,3,3)
% contourf(x2,y2,abs(z2))
hold on
plot(x3,y3,'r')
plot(x3.',y3.','b')
% xlim([-2,2])
% ylim([-2,2])


function z = slit_disc(zeta)
z = zeta;
z = (1i+z)./(1i-z);
z = z.^2;
z = 0.25*(z - 1);
% z = zeta;
% z = (1-z)./(1+z);
% z = z.^2;
% z = z - 9;
% z = sqrt(z);
% z = (1-z)./(1+z);
end