function [f,df,ddf] = smooth_transition_1_side(t,a,b,h1,h2)
tol = 1e-12;
if (b<a+tol)
    a = min(a,b) - tol;
end
f   = h1 + (h2-h1)*  fun(t,a,b);
df  =      (h2-h1)* dfun(t,a,b);
ddf =      (h2-h1)*ddfun(t,a,b);
end

function val = fun(x,a,b)
x1 = (x-a)/(b-a);
val = fun2(x1);
end

function val = dfun(x,a,b)
x1 = (x-a)/(b-a);
dx1  = 1/(b-a);
val = dfun2(x1)*dx1;
end

function val = ddfun(x,a,b)
x1 = (x-a)/(b-a);
dx1  = 1/(b-a);
val = ddfun2(x1)*dx1^2;
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