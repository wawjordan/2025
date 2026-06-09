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
% N = 11;
% xmin=0;
% xmax=1;
% dx0=0.01;
% [~,r,dx1,fh] = my_geomspace_w_refine(N,1,xmin,xmax=xmax,dx0=dx0);
% d1 = (dx0/(r-1))*log(r^(N-1));
% out = fsolve(@(x)geom_optim(x,xmin,xmax,d1*5,N),[dx0,r]);
% 
% dx02 = out(1);
% r2   = out(2);
% fh2 = @(x) xmin + dx02*((r2^(N-1)).^x - 1)/(r2-1);
% 
% figure
% hold on
% fplot(fh,[0,1],'b')
% fplot(fh2,[0,1],'r')


inputs = struct();
inputs.epsilon = 0.1;
inputs.kappa   = 0.0;
inputs.tau     = 0.0;
inputs.vinf    = 75.0;
inputs.rhoinf  = 1.0;
inputs.pinf    = 100000.0;
inputs.gamma   = 1.4;

inputs.alpha   = 5; % (degrees)
inputs.alpha   = 0; % (degrees)
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

% airfoil.plot_airfoil(true);

% [X,Y] = ndgrid(linspace(-2,2,101),linspace(0.5,4.5,101));
% [X,Y] = meshgrid(linspace(-2,2,101),linspace(-4.5,-0.5,101));
% [X,Y] = meshgrid(linspace(-5,-1,101),linspace(-2,2,101));
% u = airfoil.airfoil_x_velocity( airfoil.z_to_zeta( (X+1i*Y).*airfoil.chord + airfoil.airfoil_coords(airfoil.thetaLE) ) );
% [dudx1,dudy1] = gradient(u,0.04,0.04);
% dudx2 = airfoil.get_prim_variable_derivatives(2,1,X,Y,true);
% dudy2 = airfoil.get_prim_variable_derivatives(2,2,X,Y,true);
% 
% v = airfoil.airfoil_y_velocity( airfoil.z_to_zeta( (X+1i*Y).*airfoil.chord + airfoil.airfoil_coords(airfoil.thetaLE) ) );
% [dvdx1,dvdy1] = gradient(v,0.04,0.04);
% dvdx2 = airfoil.get_prim_variable_derivatives(3,1,X,Y,true);
% dvdy2 = airfoil.get_prim_variable_derivatives(3,2,X,Y,true);
% 
% p = airfoil.airfoil_pressure(   airfoil.z_to_zeta( (X+1i*Y).*airfoil.chord + airfoil.airfoil_coords(airfoil.thetaLE) ) );
% [dpdx1,dpdy1] = gradient(p,0.04,0.04);
% dpdx2 = airfoil.get_prim_variable_derivatives(5,1,X,Y,true);
% dpdy2 = airfoil.get_prim_variable_derivatives(5,2,X,Y,true);

% [ x, y, fx, fy ] = wake_cut_pts( airfoil_x, airfoil_y, boundary_distance, num_wake_pts, '' );

% c1 = airfoil.airfoil_curvature2(airfoil.thetaLE);
% curvature = @(t) min(c1,airfoil.airfoil_curvature2(2*pi*t));
% deriv     = @(t) abs(airfoil.airfoil_differential_arc_length(2*pi*t));
% w = @(alpha,beta,t) (1 - beta*curvature(t)).*sqrt(1+alpha^2*deriv(t).^2);
% alpha = 1;
% beta  = (1-0.5)/c1;
% phi = @(t) w(alpha,beta,t);
% N = 33;
% F1 = @(t)interp1(linspace(0,1,N),weight_grid(0,1,N,@(t)phi(t)),t,"spline");
alpha = 10; % [0,)
gamma = 0.5; % [0,1]
N = 129;

% F2 = @(d) arrayfun(@(d) TE_dist(airfoil,alpha,gamma,N,33,d),d);
% % a = [0,1,10,100,1000];
% a = (1/32) * (1/2).^(4:-1:0);
% d = F2(a);
% loglog(a,d);

F1 = airfoil_param_fun(airfoil,alpha,gamma,33);
F2 = airfoil_param_fun(airfoil,1,0,33);

hold on
fplot(F1,[0,1],'r')
fplot(F2,[0,1],'b')

F1 = airfoil_param_fun2(airfoil,alpha,gamma,33,N,0.5);
F2 = airfoil_param_fun2_new(airfoil,alpha,gamma,33,N,0.5);
F3 = airfoil_param_fun3(airfoil,alpha,gamma,33,N,0.5);

t = linspace(0,1,N);

hold on
plot(t(1:end-1),diff(F2(t)),'b')
plot(t(1:end-1),diff(F3(t)),'g')
hold off;
% [t1,Ft,~,~,~] = my_tanh_stretching_function( 0, 1, N, 1e-1, 1e-1, nan, nan, 1.0e-6 );
% fplot(phi,[0,1])
% N = 2049;
% F1 = @(t)interp1(linspace(0,1,N),weight_grid(0,1,N,@(t)phi(t)),t,"spline");
% N = 21;
% F2 = @(t)interp1(linspace(0,1,N),weight_grid(0,1,N,@(t)phi(t)),t,"spline");
% % hold on
% fplot(@(t)F1(t)-F2(t),[0,1],'r')
% fplot(@(t)F1(t),[0,1],'r')
% fplot(phi,[0,1],'b')

figure
% fplot(F1,[0,1],'r')
% hold on
fplot(F2,[0,1],'b')
hold on
fplot(F3,[0,1],'g')

airfoil.plot_airfoil(true);
[x,y] = airfoil.output_airfoil_coords1(N,F1);
plot(x,y,'r.')
[x,y] = airfoil.output_airfoil_coords1(N,F2);
plot(x,y,'b.')
[x,y] = airfoil.output_airfoil_coords1(N,F3);
plot(x,y,'g.')
xlim([-0.0579 0.0538])
hold off


%% create grid from cylinder
zeta = airfoil.z_to_zeta( (x+y*1i)*airfoil.chord + airfoil.airfoil_coords(airfoil.thetaLE) );
x_airfoil = flip(real(zeta));
y_airfoil = flip(imag(zeta));

inputs2 = struct();
inputs2.boundary_distance = 500;
inputs2.jmax              = 129;             % Number of Points Off Body
inputs2.wake_pts          = 129;             % Number of Points in Wake
inputs2.wall_spacing      = 0.0005;          % Intial spacing off the wall
inputs2.scjmax            = 0.977;           % Scaling Factor
inputs2.mu                = 0.1;             % 4th order explicit smoothing factor
inputs2.muim              = 0.5;             % Implicit smoothing factor
inputs2.alpham            = 10.0;            % Alpha scheme integration factor
inputs2.jm1               = 0.3;             % Alpha variation ramp parameter
inputs2.jm2               = 0.4;             % Alpha variation ramp parameter

% inputs2.boundary_distance = 1;
% inputs2.jmax              = 10;             % Number of Points Off Body
% inputs2.wake_pts          = 10;             % Number of Points in Wake
% inputs2.wall_spacing      = 0.01;          % Intial spacing off the wall
% inputs2.scjmax            = 1.0;           % Scaling Factor
% inputs2.mu                = 1;             % 4th order explicit smoothing factor
% inputs2.muim              = 1;             % Implicit smoothing factor
% inputs2.alpham            = 10;            % Alpha scheme integration factor
% inputs2.jm1               = 0.0;             % Alpha variation ramp parameter
% inputs2.jm2               = 0.0;             % Alpha variation ramp parameter

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

%% map it back

zeta = ( airfoil.zeta_to_z( (x2+y2*1i) ) - airfoil.airfoil_coords(airfoil.thetaLE) )/airfoil.chord;
x3 = real(zeta);
y3 = imag(zeta);
figure
hold on
plot(x3,y3,'k');
plot(x3.',y3.','k');
axis equal

hold off

function F = geom_optim(x,xmin,xmax,d1,N)
dx0  = x(1);
r    = x(2);

F(1) = (dx0/(r-1))*log(r^(N-1))       - d1;
F(2) = xmin + dx0*(r^(N-1) - 1)/(r-1) - xmax;
end

function dist = TE_dist(airfoil,alpha,gamma,N1,N2,d)
F = airfoil_param_fun(airfoil,alpha,gamma,N1);
[~,Ft,~,~,~] = my_tanh_stretching_function( 0, 1, N2, d, d, nan, nan, 1.0e-6 );
[x,y] = airfoil.output_airfoil_coords1(N2,@(t)Ft(F(t)));
dist = sqrt( (x(end)-x(end-1)).^2 + (y(end)-y(end-1)).^2 );
end

function F = airfoil_param_fun(airfoil,alpha,gamma,N)
if abs(airfoil.kappa)<1e-12
    c1 = airfoil.airfoil_curvature2(airfoil.thetaLE);
else
    [~,c1] = fminbnd(@(t)1-airfoil.airfoil_curvature2(2*pi*t),0.1,0.9);
end

curvature = @(t) min(0.99*c1,airfoil.airfoil_curvature2(2*pi*t));
deriv     = @(t) abs(airfoil.airfoil_differential_arc_length(2*pi*t));
w = @(alpha,beta,t) (1 - beta*curvature(t)).*sqrt(1+alpha^2*deriv(t).^2);
beta  = (1-gamma)/c1;
phi = @(t) w(alpha,beta,t);
F = @(t)interp1(linspace(0,1,N),weight_grid(0,1,N,@(t)phi(t)),t,"spline");
end

function F = airfoil_param_fun_pressure(airfoil,alpha,gamma,N)
c1 = airfoil.airfoil_curvature2(airfoil.thetaLE);
curvature = @(t) min(c1,airfoil.airfoil_curvature2(2*pi*t));
deriv     = @(t) abs(airfoil.airfoil_differential_arc_length(2*pi*t));
w = @(alpha,beta,t) (1 - beta*curvature(t)).*sqrt(1+alpha^2*deriv(t).^2);
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

function F = airfoil_param_fun2_new(airfoil,alpha,gamma,N1,N2,d)
dist = (1/d)/(N2-1);
% dist = 1/(N2-1);
% F1 = airfoil_param_fun(airfoil,alpha,gamma,N1);
% [Ft,~,~] = vinokur_two_sided_spacing_fcn(N2,dist,dist,true);
% F = @(t)Ft(F1(t));
[F,~,~] = vinokur_two_sided_spacing_fcn(N2,dist,dist,true);
end

function F = airfoil_param_fun3(airfoil,alpha,gamma,N1,N2,d)
dist = (1/d)/(N2-1);
off  = dist;
% dist = 1/(N2-1);
F1 = airfoil_param_fun(airfoil,alpha,gamma,N1);
% [Ft,~,~] = hermite_blend_2_vinokur(N2,dist/2,dist,off,true);
% [Ft,~,~] = hermite_blend_2_vinokur(N2,dist,dist/2,off,true);
% F = @(t)Ft(F1(t));
[F,~,~] = hermite_blend_2_vinokur(N2,dist*10,dist,off,true);
end



function F = smooth_symmetric_stretch(N,mid,spacing_end,spacing_mid,overlap)
N1 = floor(mid*N)+1;
N2 = N-N1+1;
[t1,F1,dF1,~,~] = my_tanh_stretching_function( 0, mid, N1, spacing_end, spacing_mid, nan, nan, 1.0e-6 );
[t2,F2,dF2,~,~] = my_tanh_stretching_function( mid, 1, N2, spacing_mid, spacing_end, nan, nan, 1.0e-6 );
off1 = mid*overlap;
off2 = (1-mid)*overlap;


F3 = @(t) sharp_transition(t,mid-off1,1,0) .* F1(t);
F4 = @(t) sharp_transition(t,mid+off2,0,1) .* F2(t);
if ( overlap > 1.0e-6)
    join = csape([mid-off1,mid,mid+off2],[F1(mid-off1),F1(mid),F2(mid+off2)],'complete',[dF1(mid-off1),dF2(mid+off2)]);
    F5 = @(t) sharp_interval(t,mid-off1,mid+off2,0,1) .* fnval(join,t);
    F = @(t) F3(t) + F4(t) + F5(t);
else
    F = @(t) F3(t) + F4(t);
end
end

function val = sharp_transition(x,a,v1,v2)
val = x*0 + v1;
val(x>a) = v2;
end

function val = sharp_interval(x,a,b,v1,v2)
val = x*0 + v1;
val(x>=a & x<=b) = v2;
end

% function F = smooth_symmetric_stretch(N,mid,spacing_end,spacing_mid,overlap)
% N1 = floor(mid*N)+1;
% N2 = N-N1+1;
% [~,F1,~,~] = my_tanh_stretching_function( 0, mid, N1, spacing_end, spacing_mid, nan, nan, 1.0e-6 );
% [~,F2,~,~] = my_tanh_stretching_function( mid, 1, N2, spacing_mid, spacing_end, nan, nan, 1.0e-6 );
% off1 = mid*overlap;
% off2 = (1-mid)*overlap;
% F3 = @(t) smooth_transition1(t,mid-off1,mid+off2);
% F4 = @(t) 1 - F3(t);
% F = @(t) F1(t).*F4(t) + F2(t).*F3(t);
% end

% function val = smooth_transition(x,a,b,c,d)
% % smooth transition function:
% % 1 on the closed interval [b, c] and vanishes outside the open interval (a, d)
% %   if a<b<c<d
% den1 = max(b-a,1e-12);
% den2 = max(d-c,1e-12);
% val = fun2((x-a)/den1).*fun2((d-x)/den2);
% end

function val = smooth_transition1(x,a,b)
% smooth transition function:
% 1 on the closed interval [b, c] and vanishes outside the open interval (a, d)
%   if a<b<c<d
den1 = max(b-a,1e-12);
val = fun2((x-a)/den1);
end
function val = fun1(x)
val = exp(-1./max(x,eps(1))).*(x>0);
end

function val = fun2(x)
val = fun1(x)./(fun1(x) + fun1(1-x));
end

function [t,F,dF,d,s] = my_tanh_stretching_function( xi_0, xi_1, N, d0, d1, s0, s1, tol )
%% Function: stretching_function
%
% Description: Computes a one-sided or two-sided stretching function
%
% Resource: Vinokur, M. "One One-Dimensional Stretching Functions for 
%           Finite-Difference Calculations", Journal of Computational
%           Physics, 50, 215-234, 1983.
%
% Inputs: 
%          xi_0:                   starting coordinate
%          xi_1:                   ending coordinate
%          N   :                   number of points
%          d0, d1:                 spacing at start/end of curve
%          s0, s1:                 derivative at start/end of curve
%          tol (optional):         tolerance for slopes near 1
% Outputs:
%          t:                      Parametric coordinate along a curve
%          F (optional) : function handle for the stretching function

%% Initialization

s = [s0,s1];
d = [d0,d1];

fzero_flag = true;
if (nargin == 7) % use default tolerance
    tol = 0.001;
end

L = abs(xi_1-xi_0);
xi = linspace(xi_0,xi_1,N);

if ( ~isnan(s0) ) % slope specified at start
    s0 = 1/s0;
    d0 = nan;
end

if ( ~isnan(s1) ) % slope specified at end
    s1 = 1/s1;
    d1 = nan;
end

refine_flag = false;

if ( ~isnan(d0) ) % estimate slope at start
    s0 = L/(d0*(N-1));
    refine_flag = true;
end

if ( ~isnan(d1) ) % estimate slope at end
    s1 = L/(d1*(N-1));
    refine_flag = true;
end

if ~isnan(s0) && ~isnan(s1) % two-sided stretching function
    if ( refine_flag )
        s = get_optimal_s_two_sided( xi, [s0,s1], [d0,d1], tol, fzero_flag );
        s0 = s(1);
        s1 = s(2);
    end
    A   = sqrt(s0/s1);
    B   = sqrt(s0*s1);
    [t,F,dF] = two_sided_stretching_function( xi, A, B, tol, fzero_flag );
elseif ~isnan(s0) && isnan(s1) % one-sided stretching function
    if ( refine_flag )
        s0 = get_optimal_s_one_sided( xi, s0, d0, tol, false, fzero_flag );
    end
    [t,F,dF] = one_sided_stretching_function( xi, s0, tol, false, fzero_flag );
    d(1) = F(xi(2)) - F(xi(1));
    d(2) = F(xi(end)) - F(xi(end-1));
    s(1) = s0;
    s(2) = L/(d(2)*(N-1));
elseif isnan(s0) && ~isnan(s1) % one-sided stretching function
    if ( refine_flag )
        s1 = get_optimal_s_one_sided( xi, s1, d1, tol, true, fzero_flag );
    end
    [t,F,dF] = one_sided_stretching_function( xi, s1, tol, true, fzero_flag );
    d(1) = F(xi(2)) - F(xi(1));
    d(2) = F(xi(end)) - F(xi(end-1));
    s(1) = L/(d(1)*(N-1));
    s(2) = s1;
else
    error('Error: unsupported option');
end

end
   
function [t,F,dF,d] = two_sided_stretching_function( xi, A, B, tol, fzero_flag )
%% Function: two_sided_stretching_function
%
% Description: Computes a two-sided hyperbolic tangent stretching function
%
% Resource: Vinokur, M. "One One-Dimensional Stretching Functions for 
%           Finite-Difference Calculations", Journal of Computational
%           Physics, 50, 215-234, 1983.
% Inputs: 
%          xi:       Computational coordinates along a curve
%          A:        Stretching function parameter
%          B:        Stretching function parameter
%          tol:      Tolerance to check if B is near 1
%          fzero_flag (optional): flag to refine root approximation
% Outputs:
%          t            : Parametric coordinate along a curve
%          F (optional) : function handle for the stretching function
%        

%% Compute Two-Sided Stretching Function
x1   = xi(1);
x2   = xi(end);
L    = x2-x1;
xmap = @(x) (x-x1)/L;
xmap_diff = @(x) 0*x + 1/L;
xmapinv = @(x) L*x + x1;
xmapinv_diff = @(x) 0*x + L;

if abs(B-1) < tol
    u_map = @(x) x.*(1 + 2*(B-1)*(x-0.5).*(1-x));
    u_map_diff = @(x) 6*(B-1)*x.*(1-x) - B + 2;
elseif B > 1
    if (nargin == 5)
        delta_y = hyperbolic_sine_function( B, fzero_flag );
    else
        delta_y = hyperbolic_sine_function( B );
    end
    u_map = @(x) 0.5*( 1 + tanh(delta_y*(x-0.5))/tanh(0.5*delta_y) );
    u_map_diff = @(x) 0.5*((delta_y*( 1 - tanh(delta_y*(x-0.5) ).^2 ) )/tanh(0.5*delta_y));
elseif B < 1
    if (nargin == 5)
        delta_x = sine_function( B, fzero_flag );
    else
        delta_x = sine_function( B );
    end
    u_map = @(x) 0.5*( 1 + tan(delta_x*(x-0.5))/tan(0.5*delta_x) );
    u_map_diff = @(x) 0.5*((delta_x*(1 + tan(delta_x*(x-0.5)).^2))/tan(0.5*delta_x));
end

t_map = @(x) u_map(x)./( A + (1-A)*u_map(x) );

t_map_diff = @(x) u_map_diff(x)./(A - u_map(x)*(A - 1)) + ((A-1)*u_map(x).*u_map_diff(x))./(A - (A-1)*u_map(x)).^2;

F = @(x) xmapinv( t_map( xmap( x ) ) );
dF = @(x) xmap_diff(x).*t_map_diff(xmap(x)).*xmapinv_diff(t_map(xmap(x)));
t = F(xi);

if ( nargout == 4 )
    d = [F(xi(2)) - F(xi(1)), F(xi(end)) - F(xi(end-1))];
end


end

function s = get_optimal_s_one_sided( xi, s0, d0, tol, flipped, fzero_flag )
    obj_fun = @(s) one_sided_err( xi, s, d0, tol, flipped, fzero_flag );
    s = fzero(obj_fun,s0);
end

function err = one_sided_err( xi, s0, d0, tol, flipped, fzero_flag )
    [~,~,~,d] = one_sided_stretching_function( xi, s0, tol, flipped, fzero_flag );
    err = d0 - d;
end

function s = get_optimal_s_two_sided( xi, s0, d0, tol, fzero_flag )
    obj_fun = @(s) two_sided_err( xi, s, d0, tol, fzero_flag );
    % options = optimoptions("fsolve","Display","none");
    % options = optimoptions("fsolve","FunctionTolerance",1e-30,"OptimalityTolerance",1e-30,"StepTolerance",1e-30,"Display","iter-detailed");
    options = optimoptions("fsolve","FunctionTolerance",1e-30,"OptimalityTolerance",1e-30,"StepTolerance",1e-30,"Display","none");
    s = fsolve(obj_fun,s0,options);
end

function err = two_sided_err( xi, s0, d0, tol, fzero_flag )
    % use approximate spacings from s0 if d0 is not provided
    if any(isnan(d0))
        N = length(xi);
        L = xi(N)-xi(1);
        if isnan(d0(1))
            d0(1) = L/(s0(1)*(N-1));
        end
        if isnan(d0(2))
            d0(2) = L/(s0(2)*(N-1));
        end
    end
    A   = sqrt(s0(1)/s0(2));
    B   = sqrt(s0(1)*s0(2));
    [~,~,~,d] = two_sided_stretching_function( xi, A, B, tol, fzero_flag );
    err(1) = d(1) - d0(1);
    err(2) = d(2) - d0(2);
end

function [t,F,dF,d] = one_sided_stretching_function( xi, s0, tol, flipped, fzero_flag )
%% Function: one_sided_stretching_function
%
% Description: Computes a one-sided hyperbolic tangent stretching function
%
% Resource: Vinokur, M. "One One-Dimensional Stretching Functions for 
%           Finite-Difference Calculations", Journal of Computational
%           Physics, 50, 215-234, 1983.
% Inputs: 
%          xi:       Computational coordinates along a curve
%          s0:       Stretching function parameter
%          tol:      Tolerance to check if s0 is near 1
%          flipped:  flag to show which end slope is specified on
%          fzero_flag (optional): flag to refine root approximation
% Outputs:
%          t            : Parametric coordinate along a curve
%          F (optional) : function handle for the stretching function      

%% Compute One-Sided Stretching Function
x1   = xi(1);
x2   = xi(end);
L    = x2-x1;
xmap = @(x) (x-x1)/L;
xmap_diff = @(x) 0*x + 1/L;
xmapinv = @(x) L*x + x1;
xmapinv_diff = @(x) 0*x + L;

if abs(s0-1) < tol
    t_map = @(x) x.*(1-0.5*(s0-1)*(1-x).*(2-x));
elseif s0 > 1
    if ( nargin == 4 )
        delta_y = hyperbolic_sine_function( s0, fzero_flag ) / 2;
    else
        delta_y = hyperbolic_sine_function( s0 ) / 2;
    end
    t_map = @(x) 1 + ( tanh(delta_y*(x-1))/tanh(delta_y) );
    t_map_diff = @(x) (delta_y/tanh(delta_y))*(1 - tanh(delta_y*(x-1)).^2);
elseif s0 < 1
    if ( nargin == 4 )
        delta_x = sine_function( s0, fzero_flag ) / 2;
    else
        delta_x = sine_function( s0 ) / 2;
    end
    t_map = @(x) 1 + ( tan(delta_x*(x-1))/tan(delta_x) );
    t_map_diff = @(x) (delta_x/tan(delta_x))*(1 + tan(delta_x*(x-1)).^2);
end

if ( flipped )
    F = @(x) xmapinv( 1 - t_map( 1 - xmap(x) ) );
    dF = @(x) xmap_diff(x).*t_map_diff(1-xmap(x)).*xmapinv_diff(1-t_map(1-xmap(x)));
else
    F = @(x) xmapinv( t_map( xmap(x) ) );
    dF = @(x) xmap_diff(x).*t_map_diff(xmap(x)).*xmapinv_diff(t_map(xmap(x)));
end

t = F(xi);

if ( nargout == 4 )
    if ( flipped )
        d = F(xi(end)) - F(xi(end-1));
    else
        d = F(xi(2)) - F(xi(1));
    end
end

end

function [ x ] = hyperbolic_sine_function( y,fzero_flag)
%% Function: hyperbolic_sine_function
% Description: Returns the approximate root of the function:
%                      y = sinh(x)/x
% Inputs:
%           y:   LHS of function
%           fzero_flag (optional): flag to refine root approximation
% Outputs:  
%           x:   Root of function

%% Determine Approximate Root
a = 2.7829681;
if y < a
    x = sinh_function_1(y);
else
    x = sinh_function_2(y);
end
if ( nargin == 2 )
    if ( fzero_flag )
        x = fzero(@(z)sinh(z)./z - y,x);
    end
end

end

function x = sinh_function_1(y)
y1 = y - 1;
a =  1;
b = -0.15;
c =  0.057321429;
d = -0.024907295;
e =  0.0077424461;
f = -0.0010794123;

x = sqrt(6*y1).*( a + b*y1 + c*y1.^2 + d*y1.^3 + e*y1.^4 + f*y1.^5 );
end

function x = sinh_function_2(y)
v = log(y);
w = 1./y - 0.028527431;
a = -0.02041793;
b =  0.24902722;
c =  1.9496443;
d = -2.6294547;
e =  8.56795911;

x = v + ( 1 + 1./v ).*log(2*v) + a + b*w + c*w.^2 + d*w.^3 + e*w.^4;
end

function [ x ] = sine_function( y, fzero_flag )
%% Function: sine_function
% Description: Returns the approximate root of the function:
%                      y = sin(x)/x
% Inputs:
%           y:   LHS of function
%           fzero_flag (optional): flag to refine root approximation
% Outputs:  
%           x:   Root of function

%% Determine Approximate Root
a = 0.26938972;
if y < a
    x = sine_function_1(y);
else
    x = sine_function_2(y);
end
if ( nargin == 2 )
    if ( fzero_flag )
        x = fzero(@(z)sin(z)./z - y,x);
    end
end
end

function x = sine_function_1(y)
a =  1;
b = -1;
c =  1;
d = -(1 + pi^2/6);
e =  6.794732;
f = -13.205501;
g =  11.726095;

x = pi*( a + b*y + c*y.^2 + d*y.^3 + e*y.^4 + f*y.^5 + g*y.^6 );
end


function x = sine_function_2(y)
y1 = 1 - y;
a =  1;
b =  0.15;
c =  0.057321429;
d =  0.048774238;
e = -0.053337753;
f =  0.075845134;

x = sqrt(6*y1).*( a + b*y1 + c*y1.^2 + d*y1.^3 + e*y1.^4 + f*y1.^5 );
end