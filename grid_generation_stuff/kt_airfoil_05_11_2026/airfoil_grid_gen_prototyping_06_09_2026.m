%% airfoil grid gen prototyping (05/11/2026)
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

% % airfoil.plot_airfoil(true);

alpha = 10; % [0,)
gamma = 0.5; % [0,1]
N = 129;

F1 = airfoil_param_fun(airfoil,alpha,gamma,33);
F1 = airfoil_param_fun2(airfoil,alpha,gamma,33,129,0.5);

[x_airfoil,y_airfoil] = airfoil.output_airfoil_coords1(N,F1);

inputs2 = struct();
inputs2.boundary_distance = 10;
inputs2.jmax              = 129;             % Number of Points Off Body
inputs2.wake_pts          = 129;             % Number of Points in Wake
inputs2.wall_spacing      = 0.0005;          % Intial spacing off the wall
inputs2.scjmax            = 0.977;           % Scaling Factor
inputs2.mu                = 0.1;             % 4th order explicit smoothing factor
inputs2.muim              = 0.5;             % Implicit smoothing factor
inputs2.alpham            = 1.0;            % Alpha scheme integration factor
inputs2.jm1               = 0.3;             % Alpha variation ramp parameter
inputs2.jm2               = 0.4;             % Alpha variation ramp parameter


[x2,y2] = hyperbolic_C_grid_local( x_airfoil, y_airfoil,                 ...
                                 inputs2.boundary_distance, inputs2.jmax, ...
                                 inputs2.wall_spacing, inputs2.scjmax, ...
                                 inputs2.mu, inputs2.muim, inputs2.alpham, ...
                                 inputs2.jm1, inputs2.jm2, inputs2.wake_pts );

hold on
plot(x2,y2,'k');
plot(x2.',y2.','k');
axis equal

hold off


function F = airfoil_param_fun(airfoil,alpha,gamma,N)
if abs(airfoil.kappa)<1e-12 % symmetric
    c1 = airfoil.airfoil_curvature2(airfoil.thetaLE);
else
    [~,c1] = fminbnd(@(t)1-airfoil.airfoil_curvature2(2*pi*t),0.1,0.9);
end

w_min = 1.0e-3; % minimum weight to prevent division by zero

curvature = @(t) min(c1,airfoil.airfoil_curvature2(2*pi*t));
deriv     = @(t) abs(airfoil.airfoil_differential_arc_length(2*pi*t));
w = @(alpha,beta,t) w_min + (1 - beta*curvature(t)).*sqrt(1+alpha^2*deriv(t).^2);
beta  = (1-gamma)/c1;
phi = @(t) w(alpha,beta,t);
F = @(t)interp1(linspace(0,1,N),weight_grid(0,1,N,@(t)phi(t)),t,"spline");
end

function F = airfoil_param_fun2(airfoil,alpha,gamma,N1,N2,d)
F1 = airfoil_param_fun(airfoil,alpha,gamma,N1);
dist = d/(N2-1);
[~,Ft,~,~,~] = my_tanh_stretching_function( 0, 1, N2, dist, dist, nan, nan, 1.0e-6 );
F = @(t)Ft(F1(t));
end
