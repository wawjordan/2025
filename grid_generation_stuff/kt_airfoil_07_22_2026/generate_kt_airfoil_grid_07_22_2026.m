function [ OUT, airfoil, fine_grid ] = generate_kt_airfoil_grid_07_22_2026( ...
                                                            airfoil_inputs, ...
                                                            n_body_pts,     ...
                                                            n_wake_pts,     ...
                                                            jmax,           ...
                                                            scjmax,         ...
                                                            alpham,         ...
                                                            AR_LE,          ...
                                                            AR_TE,          ...
                                                            delta_LE,       ...
                                                            delta_TE,       ...
                                                            spline_order,   ...
                                                            i_fine,         ...
                                                            j_fine )
% inputs
% jmax         - number of points off body
% n_body_pts   - number of points on the body
% n_wake_pts   - number of points in wake
% AR_LE        - target aspect ratio at the leading edge
% AR_TE        - target aspect ratio at the trailing edge
% delta_LE     - target wall spacing at leading edge
% delta_TE     - target wall spacing at trailing edge
% spline_order - order of spline for interpolating to fine grid
% i_fine       - i_max on the fine grid
% j_fine       - j_max on the fine grid

OUT = struct();
airfoil = local_airfoil_generator( airfoil_inputs );

% number of transition points for surface spacing blending
n_transition = ceil(n_body_pts/10);

imax = n_body_pts + 2*(n_wake_pts-1);
i1 = n_wake_pts;
i2 = n_wake_pts+n_body_pts-1;

wake_multiplier   = 1;
boundary_distance = 500;
mu                = 0.1; % 4th order explicit smoothing factor
muim              = 0.5; % Implicit smoothing factor

[x_airfoil,y_airfoil,theta,theta_fun,~,~,~,~,fy_wake] = ...
                                   get_airfoil_coordinates( airfoil,        ...
                                                            n_body_pts,     ...
                                                            n_transition,   ...
                                                            AR_LE,          ...
                                                            AR_TE,          ...
                                                            delta_LE,       ...
                                                            delta_TE,       ...
                                                            n_wake_pts,     ...
                                                            boundary_distance );
[x,y] = hyperbolic_C_grid_local_v2( x_airfoil,                              ...
                                    y_airfoil,                              ...
                                    n_wake_pts,                             ...
                                    jmax,                                   ...
                                    boundary_distance = boundary_distance,  ...
                                    wake_multiplier   = wake_multiplier,    ...
                                    AR_LE             = AR_LE,              ...
                                    AR_TE             = AR_TE,              ...
                                    n_extrude_layers  = 0,                  ...
                                    n_radial_passes   = 1,                  ...
                                    scjmax            = scjmax,             ...
                                    mu                = mu,                 ...
                                    muim              = muim,               ...
                                    alpham            = alpham,             ...
                                    fy_wake           = fy_wake );

OUT.airfoil_inputs              = airfoil_inputs;
OUT.base_grid                   = struct();
OUT.base_grid.n_body_pts        = n_body_pts;
OUT.base_grid.n_wake_pts        = n_wake_pts;
OUT.base_grid.jmax              = jmax;
OUT.base_grid.imax              = imax;
OUT.base_grid.i1                = i1;
OUT.base_grid.i2                = i2;
OUT.base_grid.AR_LE             = AR_LE;
OUT.base_grid.AR_TE             = AR_TE;
OUT.base_grid.delta_LE          = delta_LE;
OUT.base_grid.delta_TE          = delta_TE;
OUT.spline_order                = spline_order;
OUT.base_grid.boundary_distance = boundary_distance;
OUT.base_grid.x_airfoil         = x_airfoil;
OUT.base_grid.y_airfoil         = y_airfoil;
OUT.base_grid.x                 = x;
OUT.base_grid.y                 = y;
OUT.base_grid.theta             = theta;
OUT.base_grid.theta_fun         = theta_fun;
if (jmax>1)
    OUT.xs = spapi( {spline_order,spline_order}, ...
                    {linspace(0,1,imax),linspace(0,1,jmax)}, x );
    OUT.ys = spapi( {spline_order,spline_order}, ...
                    {linspace(0,1,imax),linspace(0,1,jmax)}, y );
end
OUT.xi_min  = 0;
OUT.xi_max  = 1;
OUT.eta_min = 0;
OUT.eta_max = 1;

if (nargout>2)
    fine_grid = struct();
    [ fine_grid.GRID, fine_grid.i1, fine_grid.i2 ] = get_fine_grid( OUT,     ...
                                                                    airfoil, ...
                                                                    i_fine,  ...
                                                                    j_fine );
end
end

function [ GRID, i1, i2 ] = get_fine_grid( OUT, airfoil, imax, jmax )
GRID      = struct();
GRID.imax = imax;
GRID.jmax = jmax;
ni = (GRID.imax-1)/(OUT.base_grid.imax-1);
nj = (GRID.jmax-1)/(OUT.base_grid.jmax-1);
i1 = (OUT.base_grid.i1-1)*ni+1;
i2 = (OUT.base_grid.i2-1)*ni+1;
% integer refinement/coarsening
if mod(ni,1)<10*eps(1) && mod(nj,1)<10*eps(1)
    xi_vec  = linspace(0,1,imax);
    eta_vec = linspace(0,1,jmax);
else
    % we need to make sure that there are points on the trailing edge
    i1       = round(i1);
    i2       = round(i2);
    off      = ([ 5, 5 ]-1)/(imax-1);
    xi1      = [ 0, i1-1, i2-1, imax-1 ];
    xi1      = xi1/xi1(end);
    xi2      = [ 0, (OUT.base_grid.i1-1)/(OUT.base_grid.imax-1), ...
                    (OUT.base_grid.i2-1)/(OUT.base_grid.imax-1), 1 ];
    [fh,~,~] = hermite_blended_piecewise_linear( xi1, xi2, off, off );
    xi_vec   = fh( linspace(0,1,imax) );
    eta_vec  =     linspace(0,1,jmax);
end
GRID.x = fnval(OUT.xs,{xi_vec,eta_vec});
GRID.y = fnval(OUT.ys,{xi_vec,eta_vec});

% enforce continuity (?)
GRID.x(1:i1,1) = 0.5*( GRID.x(1:i1,1) + GRID.x(end:-1:i2,1) );
GRID.x(end:-1:i2,1) = GRID.x(1:i1,1);

GRID.y(1:i1,1) = 0.5*( GRID.y(1:i1,1) + GRID.y(end:-1:i2,1) );
GRID.y(end:-1:i2,1) = GRID.y(1:i1,1);

tfun = OUT.base_grid.theta_fun;

% need to adjust theta slightly to make the points orthogonal
fz = @(t) airfoil.zs_from_theta(tfun(t));
options = optimset('TolFun',1e-15,'TolX',1e-17);
for i = i1+1:i2-1
    t0 = (i-i1)/(i2-i1);
    t = fzero( @(t) ( real(fz(t)) - GRID.x(i,2) )                ...
                    .* imag(airfoil.unit_normal_cmplx(tfun(t)))  ...
                  - ( imag(fz(t)) - GRID.y(i,2) )                ...
                    .* real(airfoil.unit_normal_cmplx(tfun(t))), ...
                    t0, options );
    GRID.x(i,1) = real(fz(t));
    GRID.y(i,1) = imag(fz(t));
end
end

function [ x, y, theta, theta_fun, F1, dF1, ddF1, fx_wake, fy_wake ] =       ...
                                              get_airfoil_coordinates(       ...
                                                             airfoil,        ...
                                                             n_body_pts,     ...
                                                             n_transition,   ...
                                                             AR_LE,          ...
                                                             AR_TE,          ...
                                                             delta_LE,       ...
                                                             delta_TE,       ...
                                                             n_wake_pts,     ...
                                                             boundary_distance )
% get parametric location of maximum curvature on the airfoil
% which will approximately be the leading edge
tLE = airfoil.thetaCmax/(2*pi);
% tLE = airfoil.thetaSP/(2*pi);

% total arc length
L = airfoil.airfoil_arc_length(0,1);

% get normalized arc length to point of max curvature
t0 = airfoil.airfoil_arc_length(0,tLE)/L;

% calculate spacing
uniform_spacing = 1/(n_body_pts-1);
dLE = AR_LE*delta_LE*uniform_spacing;
dTE = AR_TE*delta_TE*uniform_spacing;

% get buffer distance over which to blend the two functions
offset = ( (n_transition-1)/(n_body_pts-1) );

% generate asymmetric stretching function for airfoil surface
[F1,dF1,ddF1] = hermite_blend_2_vinokur_asym(n_body_pts,t0,dLE,dTE,offset,true);
[x,y,theta] = airfoil.output_airfoil_coords1(n_body_pts,F1);

f1 = @(theta) abs(airfoil.airfoil_differential_arc_length(theta));
f1c = chebfun(f1,[0,2*pi],'splitting','on');
s_c = cumsum(f1c)./sum(f1c,0,2*pi);
theta_s = inv(s_c,'splitting','on');
theta_fun = @(t) theta_s( F1(1-t) );
% flip to match the clockwise convention
x = flip(x);
y = flip(y);
theta = flip(theta);
if (nargin > 8 && n_wake_pts>0)
    [ ~, ~, fx_wake ] = wake_cut_pts( x, y, boundary_distance, n_wake_pts );
    fy_wake = @(t) airfoil.psi_y_from_x(0,fx_wake(t),0);
else
    fx_wake = [];
    fy_wake = [];
end

end

function airfoil = local_airfoil_generator(in)
a_ref   = sqrt(in.gamma*in.p_ref/in.rho_ref);

% nondimensionalize inputs
vinf   = in.vinf/a_ref;
rhoinf = in.rhoinf/in.rho_ref;
pinf   = in.pinf/(in.rho_ref*a_ref^2);

airfoil        = kt_airfoil( in.epsilon, in.kappa, in.tau );
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(in.alpha);
end