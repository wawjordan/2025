function [ x, y, fx, dfx, ddfx, fy, dfy, ddfy ] = wake_cut_pts( airfoil_x, airfoil_y, boundary_distance, num_wake_pts, varargin )
% defines the points in the wake along the eta-min boundary
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
addOptional(p,'fx',  [],valid_fcn_handle);
addOptional(p,'dfx', [],valid_fcn_handle);
addOptional(p,'ddfx',[],valid_fcn_handle);
addOptional(p,'fy',  [],valid_fcn_handle);
addOptional(p,'dfy', [],valid_fcn_handle);
addOptional(p,'ddfy',[],valid_fcn_handle);
parse( p, airfoil_x, airfoil_y, boundary_distance, num_wake_pts, varargin{:} );

airfoil_x         = p.Results.airfoil_x(:).';
airfoil_y         = p.Results.airfoil_y(:).';
boundary_distance = p.Results.boundary_distance;
num_wake_pts      = p.Results.num_wake_pts;
TE_loc            = p.Results.TE_loc;
TE_slope          = p.Results.TE_slope;
fx   = p.Results.fx;
dfx  = p.Results.dfx;
ddfx = p.Results.ddfx;
fy   = p.Results.fy;
dfy  = p.Results.dfy;
ddfy = p.Results.ddfy;

% 
% delta = get_delta(airfoil_x,airfoil_y);
% dt    = 1/numel(airfoil_x);
dt    = 1/num_wake_pts; %
[delta,f,df,ddf] = get_delta(airfoil_x,airfoil_y,dt);

if (isempty(fx)) % use default
    [fx,dfx,ddfx] = default_wake_x(TE_loc,delta,boundary_distance,num_wake_pts);
end

target_x_loc = TE_loc(1) + 0.01; % note: this assumes unit chord and airfoil in standard orientation

target_t = fzero(@(t)fx(t)-target_x_loc,0);
ftmp1 = @(t) zeros(size(t)) + f;
ftmp2 = @(t) zeros(size(t)) + df;
ftmp3 = @(t) zeros(size(t)) + ddf;

[fx1,dfx1,ddfx1] = hermite_blended_functions(0,target_t,ftmp1,fx,ftmp2,dfx,ftmp3,ddfx);

% also assumes small variation in y in vicinity of trailing edge
if (isempty(fy)) % use default
    [fy,dfy,ddfy] = default_wake_y(TE_loc,TE_slope,boundary_distance,fx1,dfx1,ddfx1);
end

t = linspace(0,1,num_wake_pts);

x1 = fx1(t);
y1 = fy(t);

% clf
% hold on
% plot([-t(10:-1:2),t(1:9)],diff([airfoil_x(end-9:end-1),fx(t(1:10))]),'r.')
% plot([-t(10:-1:2),t(1:9)],diff([airfoil_x(end-9:end-1),fx1(t(1:10))]),'b.')

% Piece the wake points together with the body points
x = horzcat( x1(1,end:-1:2), airfoil_x, x1(1,2:end));
y = horzcat( y1(1,end:-1:2), airfoil_y, y1(1,2:end));
end

function [delta,f,df,ddf] = get_delta(airfoil_x,airfoil_y,dt)
    % approximate spacing and 1st derivative (2nd order backward difference)
    % average deltas of upper + lower surface
    s1 = centri_param([airfoil_x(5:-1:1);airfoil_y(5:-1:1)],1);
    s2 = centri_param([airfoil_x(end-4:end);airfoil_y(end-4:end)],1);
    s = (s1+s2)/2;
    delta = s(2)-s(1);
    i = 5;
    f     = airfoil_x(1);
    df    = (  3*s(i-0) - 4*s(i-1) + 1*s(i-2) )/(2*dt);
    ddf   = (  2*s(i-0) - 5*s(i-1) + 4*s(i-2) - 1*s(i-3) )/(dt^2);
    % df  = ( 1*s(i-0) - 1*s(i-1) )/dt;
    % ddf = ( 1*s(i-0) - 2*s(i-1) + 1*s(i-2) )/(dt^2);
end 

function [fy,dfy,ddfy] = default_wake_y(TE_loc,TE_slope,boundary_distance,fx,dfx,ddfx)
% creates a linear distribution with the specified te_slope until 1/6 of the 
% wake length, a cubic fit between 1/6 and 2/3, 
% and a constant value for the last 1/3.
dy0 = TE_slope;
dy1 = 0;
body_x = TE_loc(1);
body_y = TE_loc(2);
if (abs(dy0) < sqrt(eps(1.0)))
    fy   = @(t) zeros(size(t)) + body_y;
    dfy  = @(t) zeros(size(t));
    ddfy = @(t) zeros(size(t));
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
    f1   = @(t) body_y + dy0*( fx(t) - body_x );
    df1  = @(t) dy0*dfx(t);
    ddf1 = @(t) dy0*ddfx(t);
    f2   = @(t) c(1)*fx(t).^3 + c(2)*fx(t).^2 + c(3)*fx(t) + c(4);
    df2  = @(t) ( 3*c(1)*fx(t).^2 + 2*c(2)*fx(t) + c(3) ).*dfx(t);
    ddf2 = @(t) 6*c(1)*fx(t).*(dfx(t).^2) + 3*c(1)*(fx(t).^2).*ddfx(t) ...
              + 2*c(2)*(dfx(t).^2) + 2*c(2)*fx(t).*ddfx(t) ...
              +   c(3)*ddfx(t);
    f3   = @(t) zeros(size(t)) + y1;
    df3  = @(t) zeros(size(t));
    ddf3 = @(t) zeros(size(t));
    fy   = @(t) f1(t).*(fx(t)<=x0) ...
              + f2(t).*((fx(t)>x0)&(fx(t)<=x1)) ...
              + f3(t).*(fx(t)>x1);
    dfy  = @(t) df1(t).*(fx(t)<=x0) ...
              + df2(t).*((fx(t)>x0)&(fx(t)<=x1)) ...
              + df3(t).*(fx(t)>x1);
    ddfy = @(t) ddf1(t).*(fx(t)<=x0) ...
              + ddf2(t).*((fx(t)>x0)&(fx(t)<=x1)) ...
              + ddf3(t).*(fx(t)>x1);
end
end

function [fx,dfx,ddfx] = default_wake_x(TE_loc,delta,boundary_distance,num_pts)
% Set an exponential stretching rate for the wake region
% input t ranges from 0 to 1
body_x = TE_loc(1);
[~,r,~] = my_geomspace( num_pts, TE_loc(1), xmax=boundary_distance, dx0=delta );
rN  = r^(num_pts-1);
lrN = (num_pts-1)*log(r);
dr1  = delta/(r-1);
fx   = @(t) body_x + dr1*(rN.^t - 1);
dfx  = @(t) dr1*lrN*rN.^t;
ddfx = @(t) dr1*lrN*lrN*rN.^t;
end