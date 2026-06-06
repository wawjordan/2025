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

fzero_flag = false;
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