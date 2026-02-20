function [x,varargout] = my_geomspace_w_refine(N,Nrefine,xmin,varargin)
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
% modified version to generate consistent sequences with grid refinement
% i.e. n_skip=2 means that this sequence will match the required parameters when
% every other point is removed
% Inputs:
%  xmin - starting coordinate
%  xmax - ending coordinate
%  r    - growth rate
%  N    - number of nodes
%  x    - vector of nodes
%  dx0  - initial spacing
%  dx1  - ending spacing

if (nargout>4)
    error('too many output arguments')
end
p = inputParser;
validScalarNum = @(x) isnumeric(x) && isscalar(x);
validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
addRequired(p,'N',validScalarPosInt);
addRequired(p,'Nrefine',validScalarPosInt);
addRequired(p,'xmin',validScalarNum);
addOptional(p,'xmax',nan,validScalarNum);
addOptional(p,'r',nan,validScalarPosNum);
addOptional(p,'dx0',nan,validScalarPosNum);
addOptional(p,'dx1',nan,validScalarPosNum);
parse(p,N,Nrefine,xmin,varargin{:});

% tmp_array = cell2mat(struct2cell(p.Results));
tmp_array = [p.Results.dx0,p.Results.dx1,p.Results.N,p.Results.r,p.Results.xmax,p.Results.xmin];
mask      = ~isnan(tmp_array);
if sum(mask)~=4
    error('incorrect number of input arguments')
end

dx0 = p.Results.dx0;
dx1 = p.Results.dx1;
N = p.Results.N;
r = p.Results.r;
xmax = p.Results.xmax;
xmin = p.Results.xmin;

Nrefine = p.Results.Nrefine;
if ~validScalarPosInt( (N-1)/Nrefine )
    warning('refinement factor results in non-integer number of nodes')
end


nout = max(nargout,1) - 1;
varargout = cell(nout);
if (all(mask([3,4,5,6])))
    [x,s{1},s{2},s{3}] = geomspace_1( (N-1)/Nrefine+1, xmin, xmax, r );
elseif (all(mask([1,3,4,6])))
    [x,s{1},s{2},s{3}] = geomspace_2( (N-1)/Nrefine+1, xmin, dx0, r );
elseif (all(mask([1,2,3,6])))
    [x,s{1},s{2},s{3}] = geomspace_3( (N-1)/Nrefine+1, xmin, dx0, dx1 );
elseif (all(mask([1,3,5,6])))
    [x,s{1},s{2},s{3}] = geomspace_4( (N-1)/Nrefine+1, xmin, xmax, dx0 );
else
    error('invalid combination of input arguments')
end

for k = 1:nout
    varargout{k} = s{k};
end
    
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

function [x,dx0,dx1,fh] = geomspace_1( N, xmin, xmax, r )
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
fh = @(x) eval_fun(N,xmin,xmax,dx0,r,x);
x = fh(linspace(0,1,N));
end

function [x,xmax,dx1,fh] = geomspace_2( N, xmin, dx0, r )
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
fh = @(x) eval_fun(N,xmin,xmax,dx0,r,x);
x = fh(linspace(0,1,N));
end

function [x,xmax,r,fh] = geomspace_3( N, xmin, dx0, dx1 )
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
fh = @(x) eval_fun(N,xmin,xmax,dx0,r,x);
x = fh(linspace(0,1,N));
end

function [x,r,dx1,fh] = geomspace_4( N, xmin, xmax, dx0 )
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
fh = @(x) eval_fun(N,xmin,xmax,dx0,r,x);
x = fh(linspace(0,1,N));
end