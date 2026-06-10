function [ x, y, fx, fy ] = wake_cut_pts( airfoil_x, airfoil_y, boundary_distance, num_wake_pts, varargin )
% Defines the points in the wake
p = inputParser;
validScalarNum = @(x) isnumeric(x) && isscalar(x);
validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);

validReal        = @(x) isreal(x) && ~isnan(x) && ~isinf(x);
validRealPoints  = @(x) all(arrayfun(@(x)validReal(x),x),"all") && isvector(x);
validVec         = @(x) validRealPoints(x) && numel(x)==2;
valid_fcn_handle = @(x) isa(x,'function_handle');
addRequired(p,'airfoil_x',validRealPoints);
addRequired(p,'airfoil_y',validRealPoints);
addRequired(p,'boundary_distance',validScalarPosNum);
addRequired(p,'num_wake_pts',validScalarPosInt);
addOptional(p,'TE_loc',[1;0],validVec);
addOptional(p,'TE_slope',0,validReal);
addOptional(p,'fx',[],valid_fcn_handle);
addOptional(p,'fy',[],valid_fcn_handle);
parse( p, airfoil_x, airfoil_y, boundary_distance, num_wake_pts, varargin{:} );

airfoil_x         = p.Results.airfoil_x(:).';
airfoil_y         = p.Results.airfoil_y(:).';
boundary_distance = p.Results.boundary_distance;
num_wake_pts      = p.Results.num_wake_pts;
TE_loc            = p.Results.TE_loc;
TE_slope          = p.Results.TE_slope;
fx = p.Results.fx;
fy = p.Results.fy;

delta = get_delta(airfoil_x,airfoil_y);

if (isempty(fx)||isempty(fy)) % use default
    fx = default_wake_x(TE_loc,delta,boundary_distance,num_wake_pts);
    fy = default_wake_y(TE_loc,TE_slope,boundary_distance,fx);
end


t = linspace(0,1,num_wake_pts);

x1 = fx(t);
y1 = fy(t);

% Piece the wake points together with the body points
x = horzcat( x1(1,end:-1:2), airfoil_x, x1(1,2:end));
y = horzcat( y1(1,end:-1:2), airfoil_y, y1(1,2:end));
end

function delta = get_delta(airfoil_x,airfoil_y)
    x0 = airfoil_x(1);
    x1 = airfoil_x(2);
    y0 = airfoil_y(1);
    y1 = airfoil_y(2);

    xn   = airfoil_x(end);
    xnm1 = airfoil_x(end-1);
    yn   = airfoil_y(end);
    ynm1 = airfoil_y(end-1);
    delta1 = sqrt( (x1-x0)^2 + (y1-y0)^2 );
    delta2 = sqrt( (xn-xnm1)^2 + (yn-ynm1)^2 );
    delta = (delta1 + delta2)/2;
end

function fy = default_wake_y(TE_loc,TE_slope,boundary_distance,fx)
% creates a linear distribution with the specified te_slope until 1/6 of the 
% wake length, a cubic fit between 1/6 and 2/3, 
% and a constant value for the last 1/3.
dy0 = TE_slope;
dy1 = 0;
body_x = TE_loc(1);
body_y = TE_loc(2);
if (abs(dy0) < sqrt(eps(1.0)))
    fy = @(t) 0*t + body_y; % constant
else
    x0 = 1 + (1/6)*boundary_distance;
    y0 = body_y + dy0*( x0 - body_x );
    x1 = 1 + (2/3)*boundary_distance;
    y1 = 2*y0 - body_y;
    A = [ x0^3,   x0^2, x0, 1; ...
          x1^3,   x1^2, x1, 1; ...
        3*x0^2, 2*x0,    1, 0; ...
        3*x1^2, 2*x1,    1, 0 ];
    b = [y0; y1; dy0; dy1];
    c = A\b;
    f1 = @(t) body_y + dy0*( fx(t) - body_x );
    f2 = @(t) c(1)*fx(t).^3 + c(2)*fx(t).^2 + c(3)*fx(t) + c(4);
    f3 = @(t) 0*t + y1;
    fy = @(t) f1(t).*(fx(t)<=x0) ...
            + f2(t).*((fx(t)>x0)&(fx(t)<=x1)) ...
            + f3(t).*(fx(t)>x1);
end
end

function fx = default_wake_x(TE_loc,delta,boundary_distance,num_pts)
% Set an exponential stretching rate for the wake region
% input t ranges from 0 to 1
body_x = TE_loc(1);
% alpha = epsilon_solve(TE_loc(1), boundary_distance, delta , num_pts, 1e-10, 100);
% r = 1.0 + alpha;
[~,r,~] = my_geomspace( num_pts, TE_loc(1), xmax=boundary_distance, dx0=delta );
rN = r^(num_pts-1);
fx = @(t) body_x + delta*(rN.^t - 1)/(r-1);
    % function [root] = epsilon_solve( min_val, max_val, delta, n_points, toler, n_iter)
    % % Newton-Raphson root finding to find epsilon for
    % % R^{j+1} = R^{j} + delta*(1+epsilon)^(j-2)
    % converged = false;
    % nptm2 = n_points - 2;
    % % Estimate epsilon
    % root = (max_val/delta)^(1/nptm2) - 1; 
    % for n = 1:n_iter
    %   const = (root+1)^nptm2;
    %   value = max_val - min_val - delta*(const*(root+1) - 1)/root;
    %   if ( abs(value) < toler )
    %     converged = true;
    %     break
    %   end
    %   root = root + value/(delta*( 1 + const*(root*nptm2 - 1) )/root^2);
    % end
    % if converged == false
    %   error('Root finding function has failed to converge'); 
    % end
    % end
end