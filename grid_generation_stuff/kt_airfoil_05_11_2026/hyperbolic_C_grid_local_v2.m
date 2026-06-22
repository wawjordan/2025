function [x,y] = hyperbolic_C_grid_local_v2( x_airfoil, y_airfoil, ...
                                             n_wake_pts, jmax, varargin )
p = inputParser;
validScalarNum = @(x) isnumeric(x) && isscalar(x);
validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
validScalarNonNegInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x >= 0);
validReal        = @(x) isreal(x) && ~isnan(x) && ~isinf(x);
validAlpha       = @(x) isa(x,'function_handle') || all(arrayfun(@(x)validScalarPosNum(x),x),"all");
validRealPoints  = @(x) all(arrayfun(@(x)validReal(x),x),"all") && isvector(x);
valid_fcn_handle = @(x) isa(x,'function_handle')||isempty(x);
addRequired(p,'x_airfoil',validRealPoints);
addRequired(p,'y_airfoil',validRealPoints);
addRequired(p,'n_wake_pts',validScalarNonNegInt);
addRequired(p,'jmax',validScalarPosInt);
addOptional(p,'boundary_distance',100,validScalarPosNum);
addOptional(p,'wake_multiplier',1,validScalarPosNum);
addOptional(p,'AR_LE', 1,validScalarPosNum);
addOptional(p,'AR_TE', 1,validScalarPosNum);
addOptional(p,'scjmax',0.999,validAlpha);
addOptional(p,'mu',0.1,validAlpha);
addOptional(p,'muim',0.5,validAlpha);
addOptional(p,'alpham',0.5,validAlpha);
addOptional(p,'fx_wake',[],valid_fcn_handle);
addOptional(p,'fy_wake',[],valid_fcn_handle);
parse( p, x_airfoil, y_airfoil, n_wake_pts, jmax, varargin{:} );
boundary_distance       = p.Results.boundary_distance;
wake_multiplier       = p.Results.wake_multiplier;
AR_LE       = p.Results.AR_LE;
AR_TE       = p.Results.AR_TE;
scjmax      = p.Results.scjmax;
mu          = p.Results.mu;
muim        = p.Results.muim;
alpham     = p.Results.alpham;
fx_wake     = p.Results.fx_wake;
fy_wake     = p.Results.fy_wake;

x_airfoil = x_airfoil(:).';
y_airfoil = y_airfoil(:).';

%% Derived values
% Number of points along eta_min boundary
n_body_pts = size(x_airfoil,2);
if n_wake_pts>0
    imax = 2*(n_wake_pts-1) + n_body_pts;
else
    imax = n_body_pts;
end

%% Allocate Grid Storage
x = zeros(imax,jmax);
y = zeros(imax,jmax);

if isa(alpham,'function_handle')
    [I,J] = ndgrid(1:imax,1:jmax);
    alpham = alpham(I,J);
elseif isscalar(alpham)
    alpham = alpham*ones(imax,jmax);
end

if isa(scjmax,'function_handle')
    [I,J] = ndgrid(1:imax,1:jmax);
    scjmax = scjmax(I,J);
elseif isscalar(scjmax)
    scjmax = scjmax*ones(imax,jmax);
end

[TE_slope,TE_loc] = set_TE_slope( x_airfoil, y_airfoil );

%% Create ETAMIN Boundary
[ x1, y1 ] = wake_cut_pts( x_airfoil, y_airfoil, boundary_distance, ...
                           n_wake_pts, TE_slope=TE_slope, TE_loc=TE_loc, ...
                           fx=fx_wake, fy=fy_wake );
x(:,1) = x1;
y(:,1) = y1;

%% Set Up Initial Radial Spacing
% See: Kinsey and Barth - pg 23
r_space = zeros(imax, jmax);
if n_wake_pts>0
    i1 = n_wake_pts;
    i2 = n_wake_pts+n_body_pts-1;
else
    i1 = 1;
    i2 = imax;
end
[wall_spacing,iLE] = get_wall_spacing(x(:,1),y(:,1),n_wake_pts,i1,i2,AR_LE,AR_TE,wake_multiplier);

% Stretching factor

% [~,r,~] = my_geomspace( jmax, 0, xmax=boundary_distance, dx0=wall_spacing(iLE) );
% alpha = r-1;

alpha   = zeros(imax,1);
for i = 1:imax
    [~,r,~] = my_geomspace( jmax, 0, xmax=boundary_distance, dx0=wall_spacing(i) );
    alpha(i) = r-1;
end


for j = 2:jmax
    r_space(:,j) = r_space(:,j-1)                           ...
                  + wall_spacing.*(1+alpha).^(j-2);
end




%% Set Up Scaling Factor
% Note: Sets how the grid changes from the initial spacing given by 
%       wall_spacing to cells with equal volumes
% See:  Kinsey and Barth - pg 26

scale0      = 1 - exp(log(scjmax)/(jmax-2));
% for j = 2:jmax
%     scaling(:,j) = (1-scale0)^(j-1);
% end
scale    = zeros(imax,jmax);
scale(:,1) = 1;
for j = 2:jmax
    for i = 1:imax
        % scale(i,j) = (1-scale2(i))^(j-1);
        scale(i,j) = ( 1-scale0(i,j) )^(j-1);
    end
end

vscale   = ones(imax,jmax);
% for j = 1:jmax
%     vscale(:,j) = ( 1 - (wall_spacing-min(wall_spacing))./max(wall_spacing) ).^((j-1)/(jmax-1));
% end


% directly extrude the first layer
h = r_space(:,2);
n_layers = 1;
% h = wall_spacing;
for j = 2:n_layers
    [x(:,j), y(:,j)] = extrude_surface_pts(x(:,j-1),y(:,j-1),n_wake_pts,n_wake_pts+n_body_pts,h);
    h = h*r;
end

% hold on
% plot(x(:,1:n_layers),y(:,1:n_layers),'k')
% plot(x(:,1:n_layers).',y(:,1:n_layers).','k')
% axis equal
% xlim([0.9,1.1])
%% March Grid
for n = 1:2
    fprintf('pass %d / %d\n',n,2)
    for j = 2:jmax
        [x(:,j), y(:,j)] = march_grid( jmax, x(:,j-1), y(:,j-1), ...
            mu, muim, alpham(:,j),            ...
            scale(:,j), vscale(:,j), r_space(:,j), r_space(:,j-1) );
        % if any(j ==5:40)
        %   hold on
        %   plot(x(:,1:j),y(:,1:j),'k')
        %   plot(x(:,1:j).',y(:,1:j).','k')
        %   axis equal
        %   xlim([499.9,500.1])
        %   j;
        % end
        fprintf('step %d / %d\n',j-1,jmax-1)
    end

    alpha   = zeros(imax,1);
    for i = 1:imax
        d = min(sqrt( ( x(i,end) - x(:,1) ).^2 + (y(i,end) - y(:,1) ).^2 ));
        dist = boundary_distance - d;
        if dist < 0
            sf =  d / (d+abs(dist));
        else
            sf = ( d+abs(dist) ) / d;
        end
        sf = max(min(sf,3),0.3);
        [~,r,~] = my_geomspace( jmax, 0, xmax=sf*boundary_distance, dx0=wall_spacing(i) );
        alpha(i) = r-1;
        for j = 2:jmax
            r_space(:,j) = r_space(:,j-1)                           ...
                + wall_spacing.*(1+alpha).^(j-2);
        end
    end
end
end

function [TE_slope,TE_loc] = set_TE_slope( x, y )
up_slope  = ( y(end) - y(end-1) )/( x(end) - x(end-1) );
low_slope = ( y(1)   - y(2)     )/( x(1)   - x(2)     );   
TE_slope  = ( up_slope + low_slope )/2;    
TE_loc    = 0.5*[x(1)+x(end),y(1)+y(end)];
end

function [wall_spacing,iLE] = get_wall_spacing(x,y,n_wake_pts,i1,i2,AR_LE,AR_TE,wake_multiplier,iLE)

% segment lengths
deltas = sqrt( diff(x).^2 + diff(y).^2 );

% average top and bottom
TE_spacing = 0.5*( deltas(i1) + deltas(i2-1) );

if nargin<9
    % hunt for the minimum LE spacing
    idx1 = i1;
    idx2 = i2-1;
    % search forward
    for i = i1:i2-1
        if deltas(i+1) < deltas(i)
            idx1 = i;
            break
        end
    end
    
    % search backward
    for i = i2-1:-1:i1+1
        if deltas(i-1) < deltas(i)
            idx2 = i;
            break
        end
    end

    % swap if necessary
    if (idx1>idx2)
        tmp = idx2;
        idx2 = idx1;
        idx1 = tmp;
    end

    % find minumum spacing in this range
    % iLE = idx1;
    % minval = max(deltas);
    % for i = idx1:idx2
    %     if deltas(i)<minval
    %         minval = deltas(i);
    %         iLE = i;
    %     end
    % end
    [~,iLE] = min(deltas(idx1:idx2));
    iLE = iLE + idx1;
end

% average top and bottom
LE_spacing = 0.5*(deltas(iLE+1)+deltas(iLE));

% get wall spacings
wall_spacingLE = LE_spacing/AR_LE;
wall_spacingTE = TE_spacing/AR_TE;

imax = numel(x);
wall_spacing = zeros(size(x));
wall_spacing(i1:i2) = wall_spacingTE + (wall_spacingLE-wall_spacingTE)*smooth_transition(i1:i2,i1,iLE,iLE,i2);
if ( n_wake_pts > 0)
    wall_spacing(i1-1:-1:1) = smooth_transition_1_side (1:i1-1,1,i1-1,wall_spacingTE,wall_spacingTE*wake_multiplier);
    wall_spacing(i2+1:imax) = smooth_transition_1_side (i2+1:imax,i2+1,imax,wall_spacingTE,wall_spacingTE*wake_multiplier);
end
end

function [ x_update, y_update ] = march_grid( jmax, x, y, mu, muim, alpham, scale, vscale, rj, rjm1 )
%MARCH_GRID Main driver of grid marching in the eta direction
% INPUTS
% j       : current j index             (integer)
% jmax    : maximum j index             (integer)
% x       : x coordinates at level j-1  (double size [1,imax])
% y       : y coordinates at level j-1  (double size [1,imax])
% mu      : explicit scheme dissipation (double size [1,imax])
% muim    : implicit scheme dissipation (double size [1,imax])
% alpham  : implicit scheme parameter   (double size [1,imax])
%          (alpham = 0.5 ~ trapezoid rule)
%          (alpham = 1.0 ~ backward Euler)
% scale   : volume scaling factor       (double size [1,imax])
% ds      : radial spacing              (double size [1,imax])


%% Now start building the next grid line

% Calculate the estimated grid volumes
volume = calc_volume( jmax, x, y, scale, rj, rjm1 );

% Calculate the xi metrics and estimate the eta metrics
[xxi, yxi, xeta, yeta] = grid_metrics(x, y, volume);

% Set up the matrix problem
% rhs = create_rhs( x, y, volume, xxi, yxi, xeta, yeta, alpham, mu );
rhs = create_rhs( x, y, volume, xxi, yxi, xeta, yeta, scale, mu );
lhs = create_lhs( xxi, yxi, xeta, yeta, alpham, muim );
% https://www.mathworks.com/matlabcentral/answers/46316-sparse-block-diagonal-matrix

%% Solve and update
soln = lhs\rhs;

imax = numel(x);
soln = reshape(soln,[2,imax]);
x_update = x + reshape(soln(1,:),size(x));
y_update = y + reshape(soln(2,:),size(y));

end

function [volume] = calc_volume_old( jmax, x, y, scale, rj, rjm1 )
%VOLUME Based on scaling and the local arc length, calculate the estimated
% volumes, V^{0}
imax = numel(x);
volume = zeros(size(x));
arc_length = 0;
for i = 1:imax-1
    arc_length = arc_length + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 );
end

for i = 2:imax-1
  ds = ( rj(i+1) - rjm1(i+1) + rj(i-1)- rjm1(i-1) ) / 2;
  % ds = r(i,j) - r(i,j-1);
  dl = sqrt( ( (x(i+1) - x(i-1))/2 )^2 + ( ( y(i+1) - y(i-1))/2 )^2 );
  volume(i) = ds * ( scale(i)*dl + (1-scale(i)) * arc_length / (jmax-1) );
end

i = 1;
volume(i) = volume(i+1);

i = imax;
volume(i) = volume(i-1);
end

function [volume] = calc_volume( jmax, x, y, scale, rj, rjm1 )
%VOLUME Based on scaling and the local arc length, calculate the estimated
% volumes, V^{0}
imax = numel(x);
volume = zeros(size(x));
x = x(:);
y = y(:);

dx = gradient(x(:));
dy = gradient(y(:));
mag = (rj(:) - rjm1(:))./sqrt( dx.^2 + dy.^2 );
x_ = x - dy.*mag;
y_ = y + dx.*mag;

% [x2,y2] = extrude_surface_pts(x,y,1,imax,mag.*sqrt( dx.^2 + dy.^2 ));

x1 = x(1:imax-1);  y1 = y(1:imax-1);
x2 = x(2:imax);    y2 = y(2:imax);
x3 = x_(2:imax);   y3 = y_(2:imax);
x4 = x_(1:imax-1); y4 = y_(1:imax-1);

vtmp = 0.5*abs( ( x1.*y2 + x2.*y3 + x3.*y4 + x4.*y1 ) ...
              - ( x2.*y1 + x3.*y2 + x4.*y3 + x1.*y4 ) );
volume(1) = vtmp(1);
volume(2:imax-1) = 0.5*( vtmp(1:end-1) + vtmp(2:end) );
volume(imax) = vtmp(end);
% [xxi, yxi, xeta, yeta] = grid_metrics(x, y, volume);
% volj = xxi.*yeta - yxi.*xeta;
% 
% volume = scale.*volume + (1-scale).*volj;

volume_old = calc_volume_old( jmax, x, y, scale, rj, rjm1 );
% volume = scale.*volume_old + (1-scale).*volume;
volume = volume_old;

end

function [xxi, yxi, xeta, yeta] = grid_metrics(x, y, volume)
%GRID_METRICS compute the current xi metrics and then estimate the eta
% metrics for the hyperbolic marching scheme

imax = numel(x);

xxi  = zeros(size(x));
yxi  = zeros(size(x));
xeta = zeros(size(x));
yeta = zeros(size(x));

for i = 2:imax-1            
    % xi metrics are set by the previous grid line
    xxi(i) = 0.5*( x(i+1) - x(i-1) );
    yxi(i) = 0.5*( y(i+1) - y(i-1) );

    % eta metrics are guessed from the xi metrics and volumes
    factor = xxi(i)^2 + yxi(i)^2;

    xeta(i) = -yxi(i)*volume(i)/factor;
    yeta(i) =  xxi(i)*volume(i)/factor;
end

i = 1;
xxi(i) = 0.5*( -3*x(i) + 4*x(i+1) - x(i+2) );
yxi(i) = 0.5*( -3*y(i) + 4*y(i+1) - y(i+2) );
factor = xxi(i)^2 + yxi(i)^2;
xeta(i) = -yxi(i)*volume(i)/factor;
yeta(i) =  xxi(i)*volume(i)/factor;

i = imax;
xxi(i) =-0.5*( -3*x(i) + 4*x(i-1) - x(i-2) );
yxi(i) =-0.5*( -3*y(i) + 4*y(i-1) - y(i-2) );
factor = xxi(i)^2 + yxi(i)^2;
xeta(i) = -yxi(i)*volume(i)/factor;
yeta(i) =  xxi(i)*volume(i)/factor;
        
end

function rhs = create_rhs( x, y, volume, xxi, yxi, xeta, yeta, alpha, mu )
%CREATE_RHS Creates the standard rhs matrix for volume marching plus
imax = numel(x);

if isscalar(alpha)
    alpha = alpha*ones(imax,1);
end

if isscalar(mu)
    mu = mu*ones(imax,1);
end

% See Steger & Chaussee eq's 10, 11, and 12
rhs = zeros(2*imax,1);

%% Compute RHS for Interior Points
% for i = 2:imax-1
% % Explicit fourth order dissipation
%   if i == 2
%     dissipation_x = -mu(i)*(-2*x(i-1) + 5*x(i) - 4*x(i+1) + x(i+2));
%     dissipation_y = -mu(i)*(-2*y(i-1) + 5*y(i) - 4*y(i+1) + y(i+2));
%   elseif i == imax-1
%     dissipation_x = -mu(i)*(-2*x(i+1) + 5*x(i) - 4*x(i-1) + x(i-2));
%     dissipation_y = -mu(i)*(-2*y(i+1) + 5*y(i) - 4*y(i-1) + y(i-2));
%   else
%     dissipation_x = -mu(i)*(x(i-2) - 4*x(i-1) + 6*x(i) - 4*x(i+1) + x(i+2));
%     dissipation_y = -mu(i)*(y(i-2) - 4*y(i-1) + 6*y(i) - 4*y(i+1) + y(i+2));
%   end
% 
%   xxi0  = xxi(i);
%   yxi0  = yxi(i);
%   xeta0 = xeta(i);
%   yeta0 = yeta(i);
% 
%   determ = 1/(xxi0^2+yxi0^2);
% 
%   volj = xxi0*yeta0 - yxi0*xeta0;
% 
%   vblend = alpha(i)*volume(i) + (1-alpha(i))*volj;
% 
%   rhs(2*(i-1)+1) = -yxi0*determ*vblend + dissipation_x;
%   rhs(2*(i-1)+2) =  xxi0*determ*vblend + dissipation_y;
% end

for i = 1:imax
% Explicit fourth order dissipation
  if i == 1
      dissipation_x = -mu(i)*( 3*x(i)   - 14*x(i+1) + 26*x(i+2) - 24*x(i+3) + 11*x(i+4) - 2*x(i+5) );
      dissipation_y = -mu(i)*( 3*y(i)   - 14*y(i+1) + 26*y(i+2) - 24*y(i+3) + 11*y(i+4) - 2*y(i+5) );
  elseif i == imax
      dissipation_x = -mu(i)*( 3*x(i)   - 14*x(i-1) + 26*x(i-2) - 24*x(i-3) + 11*x(i-4) - 2*x(i-5) );
      dissipation_y = -mu(i)*( 3*y(i)   - 14*y(i-1) + 26*y(i-2) - 24*y(i-3) + 11*y(i-4) - 2*y(i-5) );
  elseif i == 2
      dissipation_x = -mu(i)*( 2*x(i-1) -  9*x(i)   + 16*x(i+1) - 14*x(i+2) + 6*x(i+3)  - 1*x(i+4) );
      dissipation_y = -mu(i)*( 2*y(i-1) -  9*y(i)   + 16*y(i+1) - 14*y(i+2) + 6*y(i+3)  - 1*y(i+4) );
  elseif i == imax-1
      dissipation_x = -mu(i)*( 2*x(i+1) -  9*x(i)   + 16*x(i-1) - 14*x(i-2) + 6*x(i-3)  - 1*x(i-4) );
      dissipation_y = -mu(i)*( 2*y(i+1) -  9*y(i)   + 16*y(i-1) - 14*y(i-2) + 6*y(i-3)  - 1*y(i-4) );
  else
      dissipation_x = -mu(i)*( 1*x(i-2) -  4*x(i-1) +  6*x(i)   -  4*x(i+1) + 1*x(i+2));
      dissipation_y = -mu(i)*( 1*y(i-2) -  4*y(i-1) +  6*y(i)   -  4*y(i+1) + 1*y(i+2));
  end

  xxi0  = xxi(i);
  yxi0  = yxi(i);
  xeta0 = xeta(i);
  yeta0 = yeta(i);

  determ = 1/(xxi0^2+yxi0^2);

  volj = xxi0*yeta0 - yxi0*xeta0;

  vblend = alpha(i)*volume(i) + (1-alpha(i))*volj;

  rhs(2*(i-1)+1) = -yxi0*determ*vblend + dissipation_x;
  rhs(2*(i-1)+2) =  xxi0*determ*vblend + dissipation_y;
end



%% Compute BC
i = 1;
% rhs(2*(i-1)+1) = 0;
rhs(2*(i-1)+2) = 0;

i = imax;
% rhs(2*(i-1)+1) = 0;
rhs(2*(i-1)+2) = 0;

end

function lhs = create_lhs( xxi, yxi, xeta, yeta, alpha, muim )
%CREATE_LHS creates the LHS matrix for the standard central difference
% scheme plus fourth order implicit dissipation for smoothing
imax = numel(xxi);

if isscalar(alpha)
    alpha = alpha*ones(imax,1);
end

if isscalar(muim)
    muim = muim*ones(imax,1);
end

% n = imax;
% m = 2;
% nn = 2*imax;
% nnz = 5*(n-4)*m^2;
% entries = zeros(nnz,1);
% idx_i   = zeros(nnz,1);
% idx_j   = zeros(nnz,1);

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

  lhs(2*i-1:2*i,2*i-3:2*i-2) = -alpha(i)*BinvA/2;
  lhs(2*i-1:2*i,2*i+1:2*i+2) =  alpha(i)*BinvA/2;

end

%% Now add implicit smoothing to create pentadiagonal system
% Kinsey and Barth pg 5
% Lowest block diagonal
for i = 3:imax-1
  lhs(2*i-1,2*i-5) = lhs(2*i-1,2*i-5) + muim(i);
  lhs(2*i,  2*i-4) = lhs(2*i,  2*i-4) + muim(i);
end

% Lower block diagonal
i = 2;
lhs(2*i-1,2*i-3) = lhs(2*i-1,2*i-3) - 2*muim(i);
lhs(2*i,  2*i-2) = lhs(2*i,  2*i-2) - 2*muim(i);
for i = 3:imax-1
  lhs(2*i-1,2*i-3) = lhs(2*i-1,2*i-3) - 4*muim(i);
  lhs(2*i,  2*i-2) = lhs(2*i,  2*i-2) - 4*muim(i);
end

% Diagonal
i = 2;
lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 5*muim(i);
lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 5*muim(i);
for i = 3:imax-2
  lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 6*muim(i);
  lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 6*muim(i);
end
i = imax-1;
lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 5*muim(i);
lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 5*muim(i);

% Upper block diagonal
for i = 2:imax-2
  lhs(2*i-1,2*i+1) = lhs(2*i-1,2*i+1) - 4*muim(i);
  lhs(2*i,  2*i+2) = lhs(2*i,  2*i+2) - 4*muim(i);
end
i = imax-1;
lhs(2*i-1,2*i+1) = lhs(2*i-1,2*i+1) - 2*muim(i);
lhs(2*i,  2*i+2) = lhs(2*i,  2*i+2) - 2*muim(i);

% Last diagonal
for i = 2:imax-2
  lhs(2*i-1,2*i+3) = lhs(2*i-1,2*i+3) + muim(i);
  lhs(2*i,  2*i+4) = lhs(2*i,  2*i+4) + muim(i);
end

%% Boundary Condition
lhs(2,4) = -2;
lhs(2,6) = 1;
lhs(2*imax,2*imax-2) = -2;
lhs(2*imax,2*imax-4) = 1;

lhs = sparse(lhs);
end

function lhs = create_lhs_new( xxi, yxi, xeta, yeta, alpha, muim )
%CREATE_LHS creates the LHS matrix for the standard central difference
% scheme plus fourth order implicit dissipation for smoothing
imax = numel(xxi);

if isscalar(alpha)
    alpha = alpha*ones(imax,1);
end

if isscalar(muim)
    muim = muim*ones(imax,1);
end
muim = -muim;

% n = imax;
% m = 2;
% nn = 2*imax;
% nnz = 5*(n-4)*m^2;
% entries = zeros(nnz,1);
% idx_i   = zeros(nnz,1);
% idx_j   = zeros(nnz,1);

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

  lhs(2*i-1:2*i,2*i-3:2*i-2) = -alpha(i)*BinvA/2;
  lhs(2*i-1:2*i,2*i+1:2*i+2) =  alpha(i)*BinvA/2;

end

% Explicit fourth order dissipation
  % if i == 1
  %     dissipation_x = -mu(i)*( 3*x(i)   - 14*x(i+1) + 26*x(i+2) - 24*x(i+3) + 11*x(i+4) - 2*x(i+5) );
  %     dissipation_y = -mu(i)*( 3*y(i)   - 14*y(i+1) + 26*y(i+2) - 24*y(i+3) + 11*y(i+4) - 2*y(i+5) );
  % elseif i == imax
  %     dissipation_x = -mu(i)*( 3*x(i)   - 14*x(i-1) + 26*x(i-2) - 24*x(i-3) + 11*x(i-4) - 2*x(i-5) );
  %     dissipation_y = -mu(i)*( 3*y(i)   - 14*y(i-1) + 26*y(i-2) - 24*y(i-3) + 11*y(i-4) - 2*y(i-5) );
  % elseif i == 2
  %     dissipation_x = -mu(i)*( 2*x(i-1) -  9*x(i)   + 16*x(i+1) - 14*x(i+2) + 6*x(i+3)  - 1*x(i+4) );
  %     dissipation_y = -mu(i)*( 2*y(i-1) -  9*y(i)   + 16*y(i+1) - 14*y(i+2) + 6*y(i+3)  - 1*y(i+4) );
  % elseif i == imax-1
  %     dissipation_x = -mu(i)*( 2*x(i+1) -  9*x(i)   + 16*x(i-1) - 14*x(i-2) + 6*x(i-3)  - 1*x(i-4) );
  %     dissipation_y = -mu(i)*( 2*y(i+1) -  9*y(i)   + 16*y(i-1) - 14*y(i-2) + 6*y(i-3)  - 1*y(i-4) );
  % else
  %     dissipation_x = -mu(i)*( 1*x(i-2) -  4*x(i-1) +  6*x(i)   -  4*x(i+1) + 1*x(i+2));
  %     dissipation_y = -mu(i)*( 1*y(i-2) -  4*y(i-1) +  6*y(i)   -  4*y(i+1) + 1*y(i+2));
  % end

i = 1;
% lhs(2*i-1,2*i-1)  = lhs(2*i-1,2*i-1)  -  3*muim(i);
% lhs(2*i  ,2*i  )  = lhs(2*i  ,2*i  )  -  3*muim(i);
% lhs(2*i-1,2*i+1)  = lhs(2*i-1,2*i+1)  + 14*muim(i);
% lhs(2*i  ,2*i+2)  = lhs(2*i  ,2*i+2)  + 14*muim(i);
% lhs(2*i-1,2*i+3)  = lhs(2*i-1,2*i+3)  - 26*muim(i);
% lhs(2*i  ,2*i+4)  = lhs(2*i  ,2*i+4)  - 26*muim(i);
% lhs(2*i-1,2*i+5)  = lhs(2*i-1,2*i+5)  + 24*muim(i);
% lhs(2*i  ,2*i+6)  = lhs(2*i  ,2*i+6)  + 24*muim(i);
% lhs(2*i-1,2*i+7)  = lhs(2*i-1,2*i+7)  - 11*muim(i);
% lhs(2*i  ,2*i+8)  = lhs(2*i  ,2*i+8)  - 11*muim(i);
% lhs(2*i-1,2*i+9 ) = lhs(2*i-1,2*i+9 ) +  2*muim(i);
% lhs(2*i  ,2*i+10) = lhs(2*i  ,2*i+10) +  2*muim(i);
i = 2;
% lhs(2*i-1,2*i-3)  = lhs(2*i-1,2*i-3)  -  2*muim(i);
% lhs(2*i  ,2*i-2)  = lhs(2*i  ,2*i-2)  -  2*muim(i);
% lhs(2*i-1,2*i-1)  = lhs(2*i-1,2*i-1)  +  9*muim(i);
% lhs(2*i  ,2*i  )  = lhs(2*i  ,2*i  )  +  9*muim(i);
% lhs(2*i-1,2*i+1)  = lhs(2*i-1,2*i+1)  - 16*muim(i);
% lhs(2*i  ,2*i+2)  = lhs(2*i  ,2*i+2)  - 16*muim(i);
% lhs(2*i-1,2*i+3)  = lhs(2*i-1,2*i+3)  + 14*muim(i);
% lhs(2*i  ,2*i+4)  = lhs(2*i  ,2*i+4)  + 14*muim(i);
% lhs(2*i-1,2*i+5)  = lhs(2*i-1,2*i+5)  -  6*muim(i);
% lhs(2*i  ,2*i+6)  = lhs(2*i  ,2*i+6)  -  6*muim(i);
% lhs(2*i-1,2*i+7)  = lhs(2*i-1,2*i+7)  +  1*muim(i);
% lhs(2*i  ,2*i+8)  = lhs(2*i  ,2*i+8)  +  1*muim(i);

for i = 3:imax-2
    lhs(2*i-1,2*i-5)  = lhs(2*i-1,2*i-5)  -  1*muim(i);
    lhs(2*i  ,2*i-4)  = lhs(2*i  ,2*i-4)  -  1*muim(i);
    lhs(2*i-1,2*i-3)  = lhs(2*i-1,2*i-3)  +  4*muim(i);
    lhs(2*i  ,2*i-2)  = lhs(2*i  ,2*i-2)  +  4*muim(i);
    lhs(2*i-1,2*i-1)  = lhs(2*i-1,2*i-1)  -  6*muim(i);
    lhs(2*i  ,2*i  )  = lhs(2*i  ,2*i  )  -  6*muim(i);
    lhs(2*i-1,2*i+1)  = lhs(2*i-1,2*i+1)  +  4*muim(i);
    lhs(2*i  ,2*i+2)  = lhs(2*i  ,2*i+2)  +  4*muim(i);
    lhs(2*i-1,2*i+3)  = lhs(2*i-1,2*i+3)  -  1*muim(i);
    lhs(2*i  ,2*i+4)  = lhs(2*i  ,2*i+4)  -  1*muim(i);
end

i = imax-1;
% lhs(2*i-1,2*i-9)  = lhs(2*i-1,2*i-9)  +  1*muim(i);
% lhs(2*i  ,2*i-8)  = lhs(2*i  ,2*i-8)  +  1*muim(i);
% lhs(2*i-1,2*i-7)  = lhs(2*i-1,2*i-7)  -  6*muim(i);
% lhs(2*i  ,2*i-6)  = lhs(2*i  ,2*i-6)  -  6*muim(i);
% lhs(2*i-1,2*i-5)  = lhs(2*i-1,2*i-5)  + 14*muim(i);
% lhs(2*i  ,2*i-4)  = lhs(2*i  ,2*i-4)  + 14*muim(i);
% lhs(2*i-1,2*i-3)  = lhs(2*i-1,2*i-3)  - 16*muim(i);
% lhs(2*i  ,2*i-2)  = lhs(2*i  ,2*i-2)  - 16*muim(i);
% lhs(2*i-1,2*i-1)  = lhs(2*i-1,2*i-1)  +  9*muim(i);
% lhs(2*i  ,2*i  )  = lhs(2*i  ,2*i  )  +  9*muim(i);
% lhs(2*i-1,2*i+1)  = lhs(2*i-1,2*i+1)  -  2*muim(i);
% lhs(2*i  ,2*i+2)  = lhs(2*i  ,2*i+2)  -  2*muim(i);
% 
% 
i = imax;
% lhs(2*i-1,2*i-11) = lhs(2*i-1,2*i-11) +  2*muim(i);
% lhs(2*i  ,2*i-10) = lhs(2*i  ,2*i-10) +  2*muim(i);
% lhs(2*i-1,2*i-9)  = lhs(2*i-1,2*i-9)  - 11*muim(i);
% lhs(2*i  ,2*i-8)  = lhs(2*i  ,2*i-8)  - 11*muim(i);
% lhs(2*i-1,2*i-7)  = lhs(2*i-1,2*i-7)  + 24*muim(i);
% lhs(2*i  ,2*i-6)  = lhs(2*i  ,2*i-6)  + 24*muim(i);
% lhs(2*i-1,2*i-5)  = lhs(2*i-1,2*i-5)  - 26*muim(i);
% lhs(2*i  ,2*i-4)  = lhs(2*i  ,2*i-4)  - 26*muim(i);
% lhs(2*i-1,2*i-3)  = lhs(2*i-1,2*i-3)  + 14*muim(i);
% lhs(2*i  ,2*i-2)  = lhs(2*i  ,2*i-2)  + 14*muim(i);
% lhs(2*i-1,2*i-1)  = lhs(2*i-1,2*i-1)  -  3*muim(i);
% lhs(2*i  ,2*i  )  = lhs(2*i  ,2*i  )  -  3*muim(i);

%% Boundary Condition
lhs(2,4) = -2;
lhs(2,6) = 1;
lhs(2*imax,2*imax-2) = -2;
lhs(2*imax,2*imax-4) = 1;

lhs = sparse(lhs);
end
