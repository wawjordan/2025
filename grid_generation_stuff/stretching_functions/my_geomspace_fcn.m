function [f,out] = my_geomspace_fcn(varargin)
% geometric spacing -> each subsequent interval is r times the length of
% previous
% Inputs:
%  N    - number of nodes
%  xmin - starting coordinate
%  xmax - ending coordinate
%  r    - growth rate
%  dx0  - initial spacing
%  dx1  - ending spacing

p = inputParser;
validScalarNum    = @(x) isnumeric(x) && isscalar(x);
validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
addOptional(p,   'N',nan,validScalarPosInt);
addOptional(p,'xmin',nan,validScalarNum);
addOptional(p,'xmax',nan,validScalarNum);
addOptional(p,   'r',nan,validScalarPosNum);
addOptional(p, 'dx0',nan,validScalarPosNum);
addOptional(p, 'dx1',nan,validScalarPosNum);
parse(p,varargin{:});

N    = p.Results.N;
xmin = p.Results.xmin;
xmax = p.Results.xmax;
r    = p.Results.r;
dx0  = p.Results.dx0;
dx1  = p.Results.dx1;
% opts:
% 1: N
% 2: xmin
% 3: xmax
% 4: r
% 5: dx0
% 6: dx1

% number of possible combinations is nchoosek(6,4)=15

% these aren't very well posed or super useful, so we'll focus on the
% remaining 10

% 01) 1 2 = 3 4 5 6: [    N, xmin ] = f ( xmax,    r,  dx0,  dx1 )
% 02) 1 3 = 2 4 5 6: [    N, xmax ] = f ( xmin,    r,  dx0,  dx1 )
% 03) 1 4 = 2 3 5 6: [    N,    r ] = f ( xmin, xmax,  dx0,  dx1 )
% 04) 1 5 = 2 3 4 6: [    N,  dx0 ] = f ( xmin, xmax,    r,  dx1 )
% 05) 1 6 = 2 3 4 5: [    N,  dx1 ] = f ( xmin, xmax,    r,  dx0 )

% 06) 2 3 = 1 4 5 6: [ xmin, xmax ] = f (    N,    r,  dx0,  dx1 )
% 07) 2 4 = 1 3 5 6: [ xmin,    r ] = f (    N, xmax,  dx0,  dx1 )
% 08) 2 5 = 1 3 4 6: [ xmin,  dx0 ] = f (    N, xmax,    r,  dx1 )
% 09) 2 6 = 1 3 4 5: [ xmin,  dx1 ] = f (    N, xmax,    r,  dx0 )
% 10) 3 4 = 1 2 5 6: [ xmax,    r ] = f (    N, xmin,  dx0,  dx1 )
% 11) 3 5 = 1 2 4 6: [ xmax,  dx0 ] = f (    N, xmin,    r,  dx1 )
% 12) 3 6 = 1 2 4 5: [ xmax,  dx1 ] = f (    N, xmin,    r,  dx0 )
% 13) 4 5 = 1 2 3 6: [    r,  dx0 ] = f (    N, xmin, xmax,  dx1 )
% 14) 4 6 = 1 2 3 5: [    r,  dx1 ] = f (    N, xmin, xmax,  dx0 )
% 15) 5 6 = 1 2 3 4: [  dx0,  dx1 ] = f (    N, xmin, xmax,    r )
tmp_array = [ N, xmin, xmax, r, dx0, dx1 ];
mask      = ~isnan(tmp_array);

if sum(mask)~=4
    error('incorrect number of input arguments')
end

if( ~mask(1) )
    error('parsing for non-specified N not set up')
elseif ( all(mask([1 4 5 6])) )
    [ xmin, xmax ] = opt_06(    N,    r,  dx0,  dx1 );
elseif ( all(mask([1 3 5 6])) )
    [ xmin,    r ] = opt_07(    N, xmax,  dx0,  dx1 );
elseif ( all(mask([1 3 4 6])) )
    [ xmin,  dx0 ] = opt_08(    N, xmax,    r,  dx1 );
elseif ( all(mask([1 3 4 5])) )
    [ xmin,  dx1 ] = opt_09(    N, xmax,    r,  dx0 );
elseif ( all(mask([1 2 5 6])) )
    [ xmax,    r ] = opt_10(    N, xmin,  dx0,  dx1 );
elseif ( all(mask([1 2 4 6])) )
    [ xmax,  dx0 ] = opt_11(    N, xmin,    r,  dx1 );
elseif ( all(mask([1 2 4 5])) )
    [ xmax,  dx1 ] = opt_12(    N, xmin,    r,  dx0 );
elseif ( all(mask([1 2 3 6])) )
    [    r,  dx0 ] = opt_13(    N, xmin, xmax,  dx1 );
elseif ( all(mask([1 2 3 5])) )
    [    r,  dx1 ] = opt_14(    N, xmin, xmax,  dx0 );
elseif ( all(mask([1 2 3 4])) )
    [  dx0,  dx1 ] = opt_15(    N, xmin, xmax,    r );
end

f = @(x) eval_fun(N,xmin,xmax,dx0,r,x);

out = struct();
out.N    = N;
out.xmin = xmin;
out.xmax = xmax;
out.r    = r;
out.dx0  = dx0;
out.dx1  = dx1;
    
end

function val = eval_fun(N,xmin,xmax,dx0,r,x)
    val = arrayfun(@(x)eval_real(N,xmin,xmax,dx0,r,x),x);
    function val = eval_real(N,xmin,xmax,dx0,r,x)
        if abs(x)<10*eps(1)
            val = xmin;
        elseif abs(x-1)<N*eps(1)
            val = xmax;
        elseif abs(r-1)<N*eps(1)
            val = xmin + (N-1)*dx0*x;
        elseif (mod((N-1)*x,1)<(N-1)*eps(1))
            rNx = r.^round((N-1)*x);
            val = xmin + dx0*(rNx - 1)/(r-1);
        else
            rN = r^(N-1);
            val = xmin + dx0*(rN.^x - 1)/(r-1);
        end
    end
end

function [xmin,xmax] = opt_06( N, r, dx0, dx1 )
xmin = 0;
if (N<2)
    xmax = 0;
else
    if abs(dx0-dx1)<10*eps(1)
        xmax = (N-1)*dx0;
    else
        xmax =  dx0*(r^(N-1) - 1)/(r-1);
    end
end
end


function [xmin,r] = opt_07( N, xmax, dx0, dx1 )
if (N<2)
    xmin = xmax;
    r    = 1;
else
    if abs(dx0-dx1)<10*eps(1)
        xmin = xmax - (N-1)*dx0;
        r    = 1;
    else
        r = (dx1/dx0)^(1/(N-2));
        xmin = xmax - dx0*(r^(N-1) - 1)/(r-1);
    end
end
end

function [xmin,dx0] = opt_08( N, xmax, r, dx1 )
if (N<2)
    xmin = xmax;
    dx0  = dx1;
else
    if abs(r-1)<N*eps(1)
        dx0  = dx1;
        xmin = xmax - (N-1)*dx0;
    else
        dx0 = dx1/(r^(N-2));
        xmin = xmax - dx0*(r^(N-1) - 1)/(r-1);
    end
end
end


function [xmin,dx1] = opt_09( N, xmax, r, dx0 )
if (N<2)
    xmin = xmax;
    dx1  = dx0;
else
    if abs(r-1)<N*eps(1)
        xmin = xmax - (N-1)*dx0;
        dx1 = dx0;
    else
        xmin = xmax - dx0*(r^(N-1) - 1)/(r-1);
        dx1 = dx0*r^(N-2);
    end
end
end

function [xmax,r] = opt_10( N, xmin, dx0, dx1 )
if (N<2)
    xmax = xmin;
    r    = 1;
else
    if abs(dx0-dx1)<10*eps(1)
        xmax = xmin + (N-1)*dx0;
        r    = 1;
    else
        r = (dx1/dx0)^(1/(N-2));
        xmax = xmin + dx0*(r^(N-1) - 1)/(r-1);
    end
end
end


% 11) 3 5 = 1 2 4 6: [ xmax,  dx0 ] = f (    N, xmin,    r,  dx1 )
function [xmax,dx0] = opt_11( N, xmin, r, dx1 )
if (N<2)
    xmax = xmin;
    dx0  = dx1;
else
    if abs(r-1)<N*eps(1)
        dx0  = dx1;
        xmax = xmin + (N-1)*dx0;
    else
        dx0 = dx1/(r^(N-2));
        xmax = xmin + dx0*(r^(N-1) - 1)/(r-1);
    end
end
end



function [xmax,dx1] = opt_12( N, xmin, r, dx0 )
if (N<2)
    xmax = xmin;
    dx1  = dx0;
else
    if abs(r-1)<N*eps(1)
        xmax = xmin + (N-1)*dx0;
        dx1 = dx0;
    else
        xmax = xmin + dx0*(r^(N-1) - 1)/(r-1);
        dx1 = dx0*r^(N-2);
    end
end
end

% 13) 4 5 = 1 2 3 6: [    r,  dx0 ] = f (    N, xmin, xmax,  dx1 )
function [r,dx0] = opt_13( N, xmin, xmax, dx1 )
if (N<2)
    r    = 1;
    dx0 = dx1;
else
    delta_x = xmax - xmin;
    if (abs(delta_x/(N-1)-dx1)<10*eps(1))
        r    = 1;
        dx0 = dx1;
    else
        if delta_x/(N-1) < dx1
            r0 = 1.01;
        else
            r0 = 0.99;
        end
        fun = @(r) delta_x/dx1 - (r.^(N-1) - 1)/(r-1);
        options = optimset('FunValCheck','on');
        r = fzero(fun,r0,options);
        dx0 = dx1/(r^(N-2));
    end
end
end


function [r,dx1] = opt_14( N, xmin, xmax, dx0 )
if (N<2)
    r    = 1;
    dx1 = dx0;
else
    delta_x = xmax - xmin;
    if (abs(delta_x/(N-1)-dx0)<10*eps(1))
        r    = 1;
        dx1 = dx0;
    else
        if delta_x/(N-1) > dx0
            r0 = 1.01;
        else
            r0 = 0.99;
        end
        fun = @(r) delta_x/dx0 - (r.^(N-1) - 1)/(r-1);
        options = optimset('FunValCheck','on');
        r = fzero(fun,r0,options);
        dx1 = dx0*r^(N-2);
    end
end
end

function [dx0,dx1] = opt_15( N, xmin, xmax, r )
if (N<2)
    dx0 = 0;
    dx1 = 0;
else
    if abs(r-1)<N*eps(1)
        dx0 = (xmax - xmin)/(N-1);
        dx1 = dx0;
    else
        dx0 = (xmax - xmin)*(r-1)/(r^(N-1) - 1);
        dx1 = dx0*r^(N-2);
    end
end
end

