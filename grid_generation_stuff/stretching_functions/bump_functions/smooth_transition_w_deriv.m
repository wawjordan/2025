function [f,df,ddf] = smooth_transition_w_deriv(t,h1,h2,a,b,c,d)
% smooth transition function:
% h2 on the closed interval [b, c] and h1 outside the open interval (a, d)
%   if a<b<c<d
tol = 1e-12;
if (b<a+tol)
    a = min(a,b) - tol;
end
if (c>d+tol)
    d = max(c,d) + tol;
end
f   = h1 + (h2-h1)*  fun(t,a,b,c,d);
df  =      (h2-h1)* dfun(t,a,b,c,d);
ddf =      (h2-h1)*ddfun(t,a,b,c,d);
end

function val = fun(x,a,b,c,d)
val = fun2((x-a)/(b-a)).*fun2((d-x)/(d-c));
end

function val = dfun(x,a,b,c,d)
x1 = (x-a)/(b-a);
x2 = (d-x)/(d-c);
dx1  = 1/(b-a);
dx2 = -1/(d-c);
val = (dfun2(x1)*dx1).*fun2(x2) + fun2(x1).*(dfun2(x2)*dx2);
end

function val = ddfun(x,a,b,c,d)
x1 = (x-a)/(b-a);
x2 = (d-x)/(d-c);
dx1  = 1/(b-a);
dx2 = -1/(d-c);
val = (ddfun2(x1)*dx1^2).*fun2(x2) + (dfun2(x1)*dx1).*(dfun2(x2)*dx2) ...
    + (dfun2(x1)*dx1).*(dfun2(x2)*dx2) + fun2(x1).*(ddfun2(x2)*dx2^2);
end

function val = fun1(x)
val = zeros(size(x));
val(x>0) = exp(-1./x(x>0));
end

function val = dfun1(x)
val = zeros(size(x));
val(x>0) = fun1(x(x>0))./x(x>0).^2;
end

function val = ddfun1(x)
val = zeros(size(x));
val(x>0) = fun1(x(x>0)).*( 1 - 2*x(x>0) )./x(x>0).^4;
end

function val = fun2(x)
val = fun1(x)./(fun1(x) + fun1(1-x));
end

function val = dfun2(x)
fa   =   fun1(x);
dfa  =  dfun1(x);
fb   =   fun1(1-x);
dfb  = -dfun1(1-x);
val = ( dfa.*fb - fa.*dfb )./ ( fa + fb ).^2;
end

function val = ddfun2(x)
fa    =   fun1(x);
dfa   =  dfun1(x);
ddfa  = ddfun1(x);
fb    =   fun1(1-x);
dfb   = -dfun1(1-x);
ddfb  = ddfun1(1-x);
num   = dfa.*fb - fa.*dfb;
den   = ( fa + fb ).^2;
dnum  = ( ddfa.*fb + dfa.*dfb ) - ( dfa.*dfb + fa.*ddfb );
dden  = 2*( fa + fb ).*( dfa + dfb ); 
val   = ( dnum.*den - num.*dden )./den.^2;
end

% syms x z a b c d
% f1(z) = exp(-1./z).*heaviside(z);
% f2(z) = f1(z)./(f1(z) + f1(1-z));
% f3(x) = f2((x-a)./(b-a)).*f2((d-x)./(d-c));
% f = matlabFunction(f3);
% df = matlabFunction(diff(f3,x));
% ddf = matlabFunction(diff(f3,x,2));
% clear x z a b c d
% a = 0; b = 0.2; c = 0.8; d = 1.0;
% [f,df,ddf] = smooth_transition_w_deriv_fcn(h1,h2,a,b,c,d)
% f2 = @(t) f(t,a,b,c,d);
% df2 = @(t) df(t,a,b,c,d);
% ddf2 = @(t) ddf(t,a,b,c,d);
% figure
% hold on
% fplot(f2,[-1,2],'k')
% fplot(f1,[-1,2],'r')
% fplot(df2,[-1,2],'k')
% fplot(df1,[-1,2],'r')
% fplot(ddf2,[-1,2],'k')
% fplot(ddf1,[-1,2],'r')
% fplot(@(t)ddf1(t)-ddf2(t),[-1,2],'r')
% hold off
