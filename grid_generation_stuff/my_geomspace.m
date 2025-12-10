function [x,varargout] = my_geomspace(N,xmin,varargin)
p = inputParser;
validScalarNum = @(x) isnumeric(x) && isscalar(x);
validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
addRequired(p,'N',validScalarPosInt);
addRequired(p,'xmin',validScalarNum);
addOptional(p,'xmax',nan,validScalarNum);
addOptional(p,'r',nan,validScalarPosNum);
addOptional(p,'dx0',nan,validScalarPosNum);
addOptional(p,'dx1',nan,validScalarPosNum);
parse(p,N,xmin,varargin{:});
tmp_array = cell2mat(struct2cell(p.Results));
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


nout = max(nargout,1) - 1;
if (all(mask([3,4,5,6])))
    [x,s(1),s(2)] = geomspace_1( N, xmin, xmax, r );
    for k = 1:nout
        varargout{k} = s(k);
    end
elseif (all(mask([1,3,4,6])))
    [x,s(1),s(2)] = geomspace_2( N, xmin, dx0, r );
    for k = 1:nout
        varargout{k} = s(k);
    end
elseif (all(mask([1,2,3,4])))
    [x,s(1),s(2)] = geomspace_3( N, xmin, dx0, dx1 );
    for k = 1:nout
        varargout{k} = s(k);
    end
elseif (all(mask([1,3,5,6])))
    [x,s(1),s(2)] = geomspace_4( N, xmin, xmax, dx0 );
    for k = 1:nout
        varargout{k} = s(k);
    end
else
    error('invalid combination of input arguments')
end
    
end

function [x,dx0,dx1] = geomspace_1( N, xmin, xmax, r )
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
% Inputs:
%  xmin - starting coordinate
%  xmax - ending coordinate
%  r    - growth rate
%  N    - number of nodes
% Outputs:
%  x    - vector of nodes
%  dx0  - initial spacing
%  dx1  - ending spacing
gvec = r.^(0:N-2);
dx0 = ( xmax - xmin )/sum(gvec);
x = zeros(1,N);
x(1) = xmin;
x(N) = xmax;
for i = 2:N-1
    x(i) = x(i-1) + (r^(i-2)*dx0);
end
dx1 = r^(N-2)*dx0;
end

function [x,h,dx1] = geomspace_2( N, xmin, dx0, r )
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
x = zeros(1,N);
x(1) = xmin;
for i = 2:N
    x(i) = x(i-1) + (r^(i-2)*dx0);
end
h = x(N);
dx1 = r^(N-2)*dx0;
end

function [x,xmax,r] = geomspace_3( N, xmin, dx0, dx1 )
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
r = (dx1/dx0)^(1/(N-2));
x = zeros(1,N);
x(1) = xmin;
for i = 2:N
    x(i) = x(i-1) + (r^(i-2)*dx0);
end
xmax = x(N);
end

function [x,r,dx1] = geomspace_4( N, xmin, xmax, dx0 )
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes

S = xmax - xmin;

if S/(N-2) > dx0
    r0 = 1.01;
else
    r0 = 0.99;
end

fun = @(r) S/dx0 - sum(r.^(0:N-2));
options = optimset('FunValCheck','on');
r = fzero(fun,r0,options);

x = zeros(1,N);
x(1) = xmin;
for i = 2:N-1
    x(i) = x(i-1) + (r^(i-2)*dx0);
end
x(N) = xmax;
dx1 = x(N)-x(N-1);
end