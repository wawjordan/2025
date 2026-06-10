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