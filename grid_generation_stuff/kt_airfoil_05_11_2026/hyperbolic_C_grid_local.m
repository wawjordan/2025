function [x,y] = hyperbolic_C_grid_local( x_airfoil, y_airfoil,                 ...
                                  boundary_distance, jmax, wall_spacing,...
                                  scjmax, mu, muim, alpham, jm1, jm2,   ...
                                  wake_pts )
%% Function: hyperbolic_grid
%
% Description: Computes a grid using hyperbolic grid generation
%
% Resources: Kinsey  & Barth    1984 - Primary Reference
%            Cordova & Barth    1988 
%            Steger  & Chaussee 1980
%
% Inputs: 
%           x_airfoil:            Coordinates of airfoil
%           y_airfoil:            Coordinates of airfoil
%           boundary_distance:    Distance to outer boundary
%           jmax:                 Number of points of airfoil
%           wall_spacing:         Initial spacing of the airfoil 
%           scjmax:               Scaling Factor
%           mu:                   4th Order Explicit Smoothing Factor
%           muim:                 Implicit Smoothing Factor
%           alpham:               Alpha scheme integration factor
%           jm1:                  Alpha variation ramp parameter
%           jm2:                  Alpha variation ramp parameter
%           TE_slope:             Trailing edge slope
%           wake_pts:             Number of points in wake
%
% Outputs:
%           x:                    Coordinates of grid
%           y:                    Coordinates of grid
%

x_airfoil = x_airfoil(:).';
y_airfoil = y_airfoil(:).';

%% Derived values
% Number of points along eta_min boundary
body_pts = size(x_airfoil,2);
imax = 2*(wake_pts-1) + body_pts; 

%% Allocate Grid Storage
x = zeros(imax,jmax);
y = zeros(imax,jmax);

TE_slope = set_TE_slope( x_airfoil, y_airfoil );

%% Create ETAMIN Boundary
% [ x(:,1), y(:,1) ] = wake_cut( boundary_distance, wake_pts,     ...
%                                TE_slope, x_airfoil, y_airfoil );

[ x(:,1), y(:,1), ~, ~ ] = wake_cut_pts( x_airfoil, y_airfoil, boundary_distance, wake_pts, TE_slope=TE_slope );

%% Set Up Initial Radial Spacing
% See: Kinsey and Barth - pg 23
radial_spacing = zeros(imax, jmax);

% Stretching factor
[alpha] = epsilon( 0, boundary_distance, wall_spacing, jmax, 1e-10, ...
                   100);
[~,r,~] = my_geomspace( num_pts, TE_loc(1), xmax=boundary_distance, dx0=delta );

% FIXME: resize radial spacing to be 1,jmax
% for j = 2:jmax
%   radial_spacing(:,j) = radial_spacing(:,j-1)                           ...
%                       + wall_spacing*(1+alpha)^(j-2);
% end

sf = 2;

rad_scale = @(i) 1 + (sf-1)* ( (i-1)/(wake_pts-1) ).^4;

alpha_tmp        = ones(1,imax)*alpha;
wall_spacing_tmp = ones(1,imax)*wall_spacing;
if ( abs(sf-1) > 1e-12 )
    for i = 1:wake_pts-1
        wall_spacing_tmp(i) = rad_scale(wake_pts+1-i)*wall_spacing;
        alpha_tmp(i)        = epsilon( 0, boundary_distance, wall_spacing_tmp(i), jmax, 1e-10, 100);
    end
    for i = wake_pts+body_pts:imax
        wall_spacing_tmp(i) = rad_scale(i-(wake_pts+body_pts-2))*wall_spacing;
        alpha_tmp(i)        = epsilon( 0, boundary_distance, wall_spacing_tmp(i), jmax, 1e-10, 100);
    end
end
for j = 2:jmax
    for i = 1:imax
        radial_spacing(i,j) = radial_spacing(i,j-1)                           ...
            + wall_spacing_tmp(i)*(1+alpha_tmp(i))^(j-2);
    end
end
% [t,F,d,s] = my_tanh_stretching_function( xi_0, xi_1, N, d0, d1, s0, s1, tol )
% [rad_space,space_fun,~,~] = my_tanh_stretching_function( 0, boundary_distance, jmax, wall_spacing, (pi*boundary_distance/jmax), nan, nan, 1e-12 );
% radial_spacing = repmat(rad_space(:).',imax,1);

% Determine Delta_Eta
deta = zeros(1,jmax-1);
for j = 2:jmax
    deta(j-1) = radial_spacing(1,j) - radial_spacing(1,j-1);
end

%% Set Up Scaling Factor
% Note: Sets how the grid changes from the initial spacing given by 
%       wall_spacing to cells with equal volumes
% See:  Kinsey and Barth - pg 26

scale      = 1 - exp(log(scjmax)/(jmax-2));
scaling    = zeros(imax,jmax);
scaling(:,1) = 1;
for j = 2:jmax
    for i = 1:wake_pts-1
        scaling(i,j) = (1-scale*alpha_tmp(i)/alpha)^(j-1);
    end
    for i = wake_pts:wake_pts+body_pts-1
        scaling(i,j) = (1-scale)^(j-1);
    end
    for i = wake_pts+body_pts:imax
        scaling(i,j) = (1-scale*alpha_tmp(i)/alpha)^(j-1);
    end
end

i_1 = 1;
i_2 = wake_pts;
i_3 = wake_pts+body_pts;
i_4  = imax;
j_1  = 1;
j_2  = jmax;
s0   = 1e-2;
scaling = scaling_ij((1:imax).',(1:jmax),s0,scjmax,i_1,i_2,i_3,i_4,j_1,j_2);
wake_spacing = 0.25*wall_spacing;


i_1 = 1;
i_2 = (imax-1)/2-10;
i_3 = (imax-1)/2+2+10;
i_4  = imax;
radial_spacing = radial_spacing_ij((1:imax).',(1:jmax),wall_spacing,wake_spacing,boundary_distance,i_1,i_2,i_3,i_4,j_1,j_2);



%% March Grid
for j = 2:jmax
  [x(:,j), y(:,j)] = march_grid( mu, muim, alpham, jm1, jm2, j,         ...
                                 imax, jmax, x(:,j-1), y(:,j-1),        ...
                                 scaling, radial_spacing );
end
end

function [TE_slope] = set_TE_slope( x, y )
up_slope  = ( y(end) - y(end-1) )/( x(end) - x(end-1) );
low_slope = ( y(1)   - y(2)     )/( x(1)   - x(2)     );   
TE_slope  = ( up_slope + low_slope )/2;    

end

function val = radial_spacing_ij(i,j,wall_spacing,wake_spacing,boundary_distance,i_1,i_2,i_3,i_4,j_1,j_2)

x  = (j(:)-j_1)/(j_2-j_1);
N  = j_2;
d0 = wake_spacing - (wake_spacing-wall_spacing) ...
                       * smooth_transition(i(:),i_1,i_2,i_3,i_4);
% d0 = wake_spacing - (wake_spacing-wall_spacing) ...
%                        * (1 - smooth_transition_tanh(i(:),1e-2,i_1,i_2,i_3,i_4,0.01,0.01));
d0 = d0/boundary_distance;
val = zeros(numel(i),numel(j));
for n = 1:numel(i)
    val(n,:) = boundary_distance*basic_tanh_space(x,0,1,N,d0(n));
end
end

function val = scaling_ij(i,j,s0,scjmax,i_1,i_2,i_3,i_4,j_1,j_2)
base_scale = 1 - (1-scjmax)*basic_tanh_slope((j-j_1)/(j_2-j_1),0,1,s0);
val = 1 - (1-base_scale).*smooth_transition(i,i_1,i_2,i_3,i_4);
end

function val = basic_tanh_space(x,x0,x1,N,d0)
% [t,F,d,s] = my_tanh_stretching_function( xi_0, xi_1, N, d0, d1, s0, s1, tol )
[~,fun1,~,~] = my_tanh_stretching_function( x0, x1, N, d0, nan, nan, nan, 1e-12 );
val = fun1(x);
end

function val = basic_tanh_slope(x,x0,x1,s0)
% [t,F,d,s] = my_tanh_stretching_function( xi_0, xi_1, N, d0, d1, s0, s1, tol )
[~,fun1,~,~] = my_tanh_stretching_function( x0, x1, 2, nan, nan, s0, nan, 1e-12 );
val = fun1(x);
end

function [ x_update, y_update ] = march_grid( mu, muim, alpham, jm1, jm2, j, imax, jmax, x, y, scale, r )
%MARCH_GRID Main driver of grid marching in the eta direction

x_update = zeros(imax,1);
y_update = zeros(imax,1);

%% Calculate scaling factors for dissipation
if j <= jm1
  alpha_scale = (j-2)/(jm1-2);
  % mu_scale    = (j-2)/(jm1-2);
elseif j >= jm2
  alpha_scale = 1 - (j-jm2)/(jmax-jm2);
  mu_scale    = 1;
else
  alpha_scale = 1;
  mu_scale    = 1;
end

jm0 = jmax/6;
alpha_scale = smooth_transition(j,jm0,jm1,jm2,jmax);

alpha = 1 + (alpham - 1)*alpha_scale;

% Set scaled dissipation factors
scaled_mu   = mu_scale*mu;
scaled_muim = mu_scale*muim;

%% Now start building the next grid line
% Calculate the estimated grid volumes
[volumes] = volume(imax, jmax, j, x, y, scale, r);
% Calculate the xi metrics and estimate the eta metrics
[xxi, yxi, xeta, yeta] = grid_metrics(imax, x, y, volumes);

% Set up the matrix problem
[rhs] = create_rhs(imax, alpha, scaled_mu, volumes, xxi, yxi, xeta, yeta, x, y );

[lhs] = create_lhs(imax, alpha, scaled_muim, xxi, yxi, xeta, yeta );

%% Solve and update
soln = lhs\rhs;

for i = 1:imax
  x_update(i) = x(i) + soln(2*(i-1) + 1);
  y_update(i) = y(i) + soln(2*(i-1) + 2);
end

fprintf('step %d / %d\n',j-1,jmax-1)
end

function [ rhs ] = create_rhs( imax, alpha, mu, volume, xxi, yxi, xeta, yeta, x, y )
%CREATE_RHS Creates the standard rhs matrix for volume marching plus
% fourth order explicit smoothing

% See Steger & Chaussee eq's 10, 11, and 12
rhs = zeros(2*imax,1);

%% Compute RHS for Interior Points
for i = 2:imax-1
% Explicit fourth order dissipation
  if i == 2
    dissipation_x = -mu*(-2*x(i-1) + 5*x(i) - 4*x(i+1) + x(i+2));
    dissipation_y = -mu*(-2*y(i-1) + 5*y(i) - 4*y(i+1) + y(i+2));
  elseif i == imax-1
    dissipation_x = -mu*(-2*x(i+1) + 5*x(i) - 4*x(i-1) + x(i-2));
    dissipation_y = -mu*(-2*y(i+1) + 5*y(i) - 4*y(i-1) + y(i-2));
  else
    dissipation_x = -mu*(x(i-2) - 4*x(i-1) + 6*x(i) - 4*x(i+1) + x(i+2));
    dissipation_y = -mu*(y(i-2) - 4*y(i-1) + 6*y(i) - 4*y(i+1) + y(i+2));
  end

  xxi0  = xxi(i,1);
  yxi0  = yxi(i,1);
  xeta0 = xeta(i,1);
  yeta0 = yeta(i,1);

  determ = 1/(xxi0^2+yxi0^2);

  volj = xxi0*yeta0 - yxi0*xeta0;

% FIXME: there may be a parenthesis mismatch here
  rhs(2*(i-1)+1,1) = -yxi0*determ*(alpha*volume(i) + (1-alpha)*volj) ...
                   + dissipation_x;
  rhs(2*(i-1)+2,1) =  xxi0*determ*(alpha*volume(i) + (1-alpha)*volj) ...
                   + dissipation_y;
end

%% Compute BC
i = 1;
rhs(2*(i-1)+1,1) = 0;
rhs(2*(i-1)+2,1) = 0;

i = imax;
rhs(2*(i-1)+1,1) = 0;
rhs(2*(i-1)+2,1) = 0;

end

function [ lhs ] = create_lhs( imax, alpha, muim, xxi, yxi, xeta, yeta )
%CREATE_LHS creates the LHS matrix for the standard central difference
% scheme plus fourth order implicit dissipation for smoothing

% Create square matrix with 1 on the diagonal
lhs = eye(2*imax);

%% Set LHS Based on Grid Type
%% Apply Steger-Chaussee hyperbolic grid marching, eq 12
% Kinsey and Barth pg 5
for i = 2:imax-1
  xxi0 = xxi(i);
  yxi0 = yxi(i);
  jacobian = (xxi0^2+yxi0^2);
  xeta0 = xeta(i)/jacobian;
  yeta0 = yeta(i)/jacobian;

  BinvA = [ xxi0*xeta0 - yxi0*yeta0, xxi0*yeta0 + yxi0*xeta0; ...
            xeta0*yxi0 + yeta0*xxi0, yxi0*yeta0 - xeta0*xxi0 ];

  lhs(2*i-1:2*i,2*i-3:2*i-2) = -alpha*BinvA/2;
  lhs(2*i-1:2*i,2*i+1:2*i+2) =  alpha*BinvA/2;

end

%% Now add implicit smoothing to create pentadiagonal system
% Kinsey and Barth pg 5
% Lowest block diagonal
for i = 3:imax-1
  lhs(2*i-1,2*i-5) = lhs(2*i-1,2*i-5) + muim;
  lhs(2*i,  2*i-4) = lhs(2*i,  2*i-4) + muim;
end

% Lower block diagonal
i = 2;
lhs(2*i-1,2*i-3) = lhs(2*i-1,2*i-3) - 2*muim;
lhs(2*i,  2*i-2) = lhs(2*i,  2*i-2) - 2*muim;
for i = 3:imax-1
  lhs(2*i-1,2*i-3) = lhs(2*i-1,2*i-3) - 4*muim;
  lhs(2*i,  2*i-2) = lhs(2*i,  2*i-2) - 4*muim;
end

% Diagonal
i = 2;
lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 5*muim;
lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 5*muim;
for i = 3:imax-2
  lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 6*muim;
  lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 6*muim;
end
i = imax-1;
lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 5*muim;
lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 5*muim;

% Upper block diagonal
for i = 2:imax-2
  lhs(2*i-1,2*i+1) = lhs(2*i-1,2*i+1) - 4*muim;
  lhs(2*i,  2*i+2) = lhs(2*i,  2*i+2) - 4*muim;
end
i = imax-1;
lhs(2*i-1,2*i+1) = lhs(2*i-1,2*i+1) - 2*muim;
lhs(2*i,  2*i+2) = lhs(2*i,  2*i+2) - 2*muim;

% Last diagonal
for i = 2:imax-2
  lhs(2*i-1,2*i+3) = lhs(2*i-1,2*i+3) + muim;
  lhs(2*i,  2*i+4) = lhs(2*i,  2*i+4) + muim;
end

%% Boundary Condition
lhs(2,4) = -2;
lhs(2,6) = 1;
lhs(2*imax,2*imax-2) = -2;
lhs(2*imax,2*imax-4) = 1;

end

function [volume] = volume( imax, jmax, j, x, y, scale, r )
%VOLUME Based on scaling and the local arc length, calculate the estimated
% volumes, V^{0}

volume = zeros(1,imax);
arc_length = 0;
for i = 1:imax-1
    arc_length = arc_length + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 );
end

for i = 2:imax-1
  ds = ( r(i+1,j) - r(i+1,j-1) + r(i-1,j)- r(i-1,j-1) ) / 2;
  % ds = r(i,j) - r(i,j-1);
  dl = sqrt( ( (x(i+1) - x(i-1))/2 )^2 + ( ( y(i+1) - y(i-1))/2 )^2 );
  volume(1,i) = ds * ( scale(i,j)*dl + (1-scale(i,j)) * arc_length / (jmax-1) );
end

i = 1;             
volume(1,i) = volume(1,i+1);

i = imax;             
volume(1,i) = volume(1,i-1);
end

function [xxi, yxi, xeta, yeta] = grid_metrics(imax, x, y, volumes)
%GRID_METRICS compute the current xi metrics and then estimate the eta
% metrics for the hyperbolic marching scheme

xxi  = zeros(imax,1);
yxi  = zeros(imax,1);
xeta = zeros(imax,1);
yeta = zeros(imax,1);
for i = 2:imax-1            
    % xi metrics are set by the previous grid line
    xxi(i,1) = 0.5*( x(i+1,1) - x(i-1,1) );
    yxi(i,1) = 0.5*( y(i+1,1) - y(i-1,1) );

    % eta metrics are guessed from the xi metrics and volumes
    factor = xxi(i,1)^2 + yxi(i,1)^2;

    xeta(i,1) = -yxi(i,1)*volumes(1,i)/factor;
    yeta(i,1) =  xxi(i,1)*volumes(1,i)/factor;
end     
        
end

function val = smooth_transition(x,a,b,c,d)
% smooth transition function:
% 1 on the closed interval [b, c] and vanishes outside the open interval (a, d)
%   if a<b<c<d
den1 = max(b-a,1e-12);
den2 = max(d-c,1e-12);
val = fun2((x-a)/den1).*fun2((d-x)/den2);
end

function val = fun1(x)
val = exp(-1./x).*(x>0);
end

function val = fun2(x)
val = fun1(x)./(fun1(x) + fun1(1-x));
end

function [t,F,d,s] = my_tanh_stretching_function( xi_0, xi_1, N, d0, d1, s0, s1, tol )
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
    [t,F] = two_sided_stretching_function( xi, A, B, tol, fzero_flag );
elseif ~isnan(s0) && isnan(s1) % one-sided stretching function
    if ( refine_flag )
        s0 = get_optimal_s_one_sided( xi, s0, d0, tol, false, fzero_flag );
    end
    [t,F] = one_sided_stretching_function( xi, s0, tol, false, fzero_flag );
    d(1) = F(xi(2)) - F(xi(1));
    d(2) = F(xi(end)) - F(xi(end-1));
    s(1) = s0;
    s(2) = L/(d(2)*(N-1));
elseif isnan(s0) && ~isnan(s1) % one-sided stretching function
    if ( refine_flag )
        s1 = get_optimal_s_one_sided( xi, s1, d1, tol, true, fzero_flag );
    end
    [t,F] = one_sided_stretching_function( xi, s1, tol, true, fzero_flag );
    d(1) = F(xi(2)) - F(xi(1));
    d(2) = F(xi(end)) - F(xi(end-1));
    s(1) = L/(d(1)*(N-1));
    s(2) = s1;
else
    error('Error: unsupported option');
end

end
   
function [t,F,d] = two_sided_stretching_function( xi, A, B, tol, fzero_flag )
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
xmapinv = @(x) L*x + x1;

if abs(B-1) < tol
    u_map = @(x) x.*(1 + 2*(B-1)*(x-0.5).*(1-x));
elseif B > 1
    if (nargin == 5)
        delta_y = hyperbolic_sine_function( B, fzero_flag );
    else
        delta_y = hyperbolic_sine_function( B );
    end
    u_map = @(x) 0.5*( 1 + tanh(delta_y*(x-0.5))/tanh(0.5*delta_y) );
elseif B < 1
    if (nargin == 5)
        delta_x = sine_function( B, fzero_flag );
    else
        delta_x = sine_function( B );
    end
    u_map = @(x) 0.5*( 1 + tan(delta_x*(x-0.5))/tan(0.5*delta_x) );
end

t_map = @(x) u_map(x)./( A + (1-A)*u_map(x) );

F = @(x) xmapinv( t_map( xmap( x ) ) );
t = F(xi);

if ( nargout == 3 )
    d = [F(xi(2)) - F(xi(1)), F(xi(end)) - F(xi(end-1))];
end


end

function s = get_optimal_s_one_sided( xi, s0, d0, tol, flipped, fzero_flag )
    obj_fun = @(s) one_sided_err( xi, s, d0, tol, flipped, fzero_flag );
    s = fzero(obj_fun,s0);
end

function err = one_sided_err( xi, s0, d0, tol, flipped, fzero_flag )
    [~,~,d] = one_sided_stretching_function( xi, s0, tol, flipped, fzero_flag );
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
    [~,~,d] = two_sided_stretching_function( xi, A, B, tol, fzero_flag );
    err(1) = d(1) - d0(1);
    err(2) = d(2) - d0(2);
end

function [t,F,d] = one_sided_stretching_function( xi, s0, tol, flipped, fzero_flag )
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
xmapinv = @(x) L*x + x1;

if abs(s0-1) < tol
    t_map = @(x) x.*(1-0.5*(s0-1)*(1-x).*(2-x));
elseif s0 > 1
    if ( nargin == 4 )
        delta_y = hyperbolic_sine_function( s0, fzero_flag ) / 2;
    else
        delta_y = hyperbolic_sine_function( s0 ) / 2;
    end
    t_map = @(x) 1 + ( tanh(delta_y*(x-1))/tanh(delta_y) );
elseif s0 < 1
    if ( nargin == 4 )
        delta_x = sine_function( s0, fzero_flag ) / 2;
    else
        delta_x = sine_function( s0 ) / 2;
    end
    t_map = @(x) 1 + ( tan(delta_x*(x-1))/tan(delta_x) );
end

if ( flipped )
    F = @(x) xmapinv( 1 - t_map( 1 - xmap(x) ) );
else
    F = @(x) xmapinv( t_map( xmap(x) ) );
end

t = F(xi);

if ( nargout == 3 )
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
