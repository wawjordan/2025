function [ x_update, y_update ] = march_grid_fast( imax, jmax, x, y, mu, muim, alpham, scale, rj, rjm1 )
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

volume = calc_volume( imax, jmax, x, y, scale, rj, rjm1 );

% Calculate the xi metrics and estimate the eta metrics
[xxi, yxi, xeta, yeta] = grid_metrics( imax, x, y, volume);

% Set up the matrix problem
rhs = create_rhs( imax, x, y, volume, xxi, yxi, xeta, yeta, scale, mu );
lhs = create_lhs( imax, xxi, yxi, xeta, yeta, alpham, muim );

%% Solve and update
soln = lhs\rhs;

soln = reshape(soln,[2,imax]);
x_update = x + reshape(soln(1,:),[imax,1]);
y_update = y + reshape(soln(2,:),[imax,1]);

end

function [volume] = calc_volume( imax, jmax, x, y, scale, rj, rjm1 )
%VOLUME Based on scaling and the local arc length, calculate the estimated
% volumes, V^{0}
volume = zeros(imax,1);
arc_length = 0;
for i = 1:imax-1
    arc_length = arc_length + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 );
end

for i = 2:imax-1
  ds = ( rj(i+1) - rjm1(i+1) + rj(i-1)- rjm1(i-1) ) / 2;
  dl = sqrt( ( (x(i+1) - x(i-1))/2 )^2 + ( ( y(i+1) - y(i-1))/2 )^2 );
  volume(i) = ds * ( scale(i)*dl + (1-scale(i)) * arc_length / (jmax-1) );
end

i = 1;
volume(i) = volume(i+1);

i = imax;
volume(i) = volume(i-1);
end

function [xxi, yxi, xeta, yeta] = grid_metrics( imax, x, y, volume)
%GRID_METRICS compute the current xi metrics and then estimate the eta
% metrics for the hyperbolic marching scheme

xxi  = zeros(imax,1);
yxi  = zeros(imax,1);
xeta = zeros(imax,1);
yeta = zeros(imax,1);

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

function rhs = create_rhs( imax, x, y, volume, xxi, yxi, xeta, yeta, alpha, mu )
%CREATE_RHS Creates the standard rhs matrix for volume marching plus

% See Steger & Chaussee eq's 10, 11, and 12
rhs = zeros(2*imax,1);

%% Compute RHS for Interior Points
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
rhs(2*(i-1)+1) = 0;
rhs(2*(i-1)+2) = 0;

i = imax;
rhs(2*(i-1)+1) = 0;
rhs(2*(i-1)+2) = 0;

end

function lhs = create_lhs( imax, xxi, yxi, xeta, yeta, alpha, muim )
%CREATE_LHS creates the LHS matrix for the standard central difference
% scheme plus fourth order implicit dissipation for smoothing

% Create square matrix with 1 on the diagonal
lhs_tmp = eye(2*imax);

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

  lhs_tmp(2*i-1:2*i,2*i-3:2*i-2) = -alpha(i)*BinvA/2;
  lhs_tmp(2*i-1:2*i,2*i+1:2*i+2) =  alpha(i)*BinvA/2;

end



%% Now add implicit smoothing to create pentadiagonal system
% Kinsey and Barth pg 5
% Lowest block diagonal
for i = 3:imax-1
  lhs_tmp(2*i-1,2*i-5) = lhs_tmp(2*i-1,2*i-5) + muim(i);
  lhs_tmp(2*i,  2*i-4) = lhs_tmp(2*i,  2*i-4) + muim(i);
end

% Lower block diagonal
i = 2;
lhs_tmp(2*i-1,2*i-3) = lhs_tmp(2*i-1,2*i-3) - 2*muim(i);
lhs_tmp(2*i,  2*i-2) = lhs_tmp(2*i,  2*i-2) - 2*muim(i);
for i = 3:imax-1
  lhs_tmp(2*i-1,2*i-3) = lhs_tmp(2*i-1,2*i-3) - 4*muim(i);
  lhs_tmp(2*i,  2*i-2) = lhs_tmp(2*i,  2*i-2) - 4*muim(i);
end

% Diagonal
i = 2;
lhs_tmp(2*i-1,2*i-1) = lhs_tmp(2*i-1,2*i-1) + 5*muim(i);
lhs_tmp(2*i,  2*i)   = lhs_tmp(2*i,  2*i)   + 5*muim(i);
for i = 3:imax-2
  lhs_tmp(2*i-1,2*i-1) = lhs_tmp(2*i-1,2*i-1) + 6*muim(i);
  lhs_tmp(2*i,  2*i)   = lhs_tmp(2*i,  2*i)   + 6*muim(i);
end
i = imax-1;
lhs_tmp(2*i-1,2*i-1) = lhs_tmp(2*i-1,2*i-1) + 5*muim(i);
lhs_tmp(2*i,  2*i)   = lhs_tmp(2*i,  2*i)   + 5*muim(i);

% Upper block diagonal
for i = 2:imax-2
  lhs_tmp(2*i-1,2*i+1) = lhs_tmp(2*i-1,2*i+1) - 4*muim(i);
  lhs_tmp(2*i,  2*i+2) = lhs_tmp(2*i,  2*i+2) - 4*muim(i);
end
i = imax-1;
lhs_tmp(2*i-1,2*i+1) = lhs_tmp(2*i-1,2*i+1) - 2*muim(i);
lhs_tmp(2*i,  2*i+2) = lhs_tmp(2*i,  2*i+2) - 2*muim(i);

% Last diagonal
for i = 2:imax-2
  lhs_tmp(2*i-1,2*i+3) = lhs_tmp(2*i-1,2*i+3) + muim(i);
  lhs_tmp(2*i,  2*i+4) = lhs_tmp(2*i,  2*i+4) + muim(i);
end

%% Boundary Condition
lhs_tmp(2,4) = -2;
lhs_tmp(2,6) = 1;
lhs_tmp(2*imax,2*imax-2) = -2;
lhs_tmp(2*imax,2*imax-4) = 1;

lhs = sparse(lhs_tmp);
end