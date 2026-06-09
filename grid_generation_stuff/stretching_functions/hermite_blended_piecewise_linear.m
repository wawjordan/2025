function [f,df,ddf] = hermite_blended_piecewise_linear(x,y,delta_minus,delta_plus)
% piecewise linear with quintic hermite smoothing sections
%   N : number of nodes
% x,y : (N) data pairs of a piecewise linear function
% delta_minus/plus : (N-2) offset around nodes where the quintic will be used
%
% x1          (-) x2  (+)     (-) xN-1 (+)           xN
%  |-----------+---|---+- ... -+---|----+------------|
%               m1     p1      mN-2     pN-2
%
% this yields (N-1) + (N-2) = 2N-3 intervals

N = numel(x);

% compute segment slopes
si = diff(y(:))./diff(x(:));

%


% generate the knot sequence and evaluate values at offsets
xi = zeros(2*N-2,1);
yi = zeros(2*N-2,1);
xi(1) = x(1);
yi(1) = y(1);
cnt = 1;
for i = 2:N-1
    cnt = cnt + 1;
    xi(cnt) = x(i)-delta_minus(i-1);
    yi(cnt) = y(i) - si(i-1)*delta_minus(i-1);
    cnt = cnt + 1;
    xi(cnt) = x(i)+delta_plus(i-1);
    yi(cnt) = y(i) + si(i)  *delta_plus(i-1);
end
xi(end) = x(end);
yi(end) = y(end);

% fval = get_fval(0.1,xi,yi,si);
f   = @(t) arrayfun(@(t)fun(  t,xi,yi,si),t);
df  = @(t) arrayfun(@(t)dfun( t,xi,yi,si),t);
ddf = @(t) arrayfun(@(t)ddfun(t,xi,yi,si),t);

% compute values at minus/plus locations


end

function fval = fun(t,xi,yi,si)
[fval,~,~] = get_fval(t,xi,yi,si);
end
function dfval = dfun(t,xi,yi,si)
[~,dfval,~] = get_fval(t,xi,yi,si);
end
function ddfval = ddfun(t,xi,yi,si)
[~,~,ddfval] = get_fval(t,xi,yi,si);
end

function [fval,dfval,ddfval] = get_fval(t,xi,yi,si)

si = [si(:).';si(:).'];
si = si(:);

if ( t<=xi(1) )
    fval   = yi(1);
    dfval  = 0;
    ddfval = 0;
    return
elseif(t>=xi(end))
    fval   = yi(end);
    dfval  = 0;
    ddfval = 0;
    return
end

% iterate through the list
% linear segments are odd, quintic segments are even
n_bins = numel(xi);
i = 2;
idx = 0;
found = false;
while i<=n_bins
    if (t>xi(i))&&(t<xi(i+1))
        idx = i;
        found = true;
    end
    if (found)
        break;
    end
    i = i + 2;
end

if ( found )
    x1 = xi(idx);
    x2 = xi(idx+1);
    y1 = yi(idx);
    y2 = yi(idx+1);
    s1 = si(i);
    s2 = si(i+2);
    % s1 = si((idx-1)/2+1);
    % s2 = si((idx+1)/2);
    [fval,dfval,ddfval] = hermite_blend_segments(t,x1,x2,y1,y2,s1,s2);
else
    idx = 0;
    found = false;
    i = 1;
    while i<=n_bins
        if (t>=xi(i))&&(t<=xi(i+1))
            idx = i;
            found = true;
        end
        if (found)
            break;
        end
        i = i + 2;
    end
    if (~found)
        error('bin not found')
    end
    fval = yi(idx) + si(i)*(t-xi(idx));
    dfval = zeros(size(t)) + si(i);
    ddfval = zeros(size(t));
end


end