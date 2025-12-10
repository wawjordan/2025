function x = eriksson_function(xmin,xmax,N,x0,alpha)
xi = linspace(0,1,N);
cond = ( (xi>=0)&(xi<x0(1)) );
d1 = (x0(1) - 0);
d2 = (xi - 0);

x = d1*f(d2/d1,alpha).*cond;

Nlines = length(x0);

for i = 2:Nlines
    cond = ( ( xi>=x0(i-1) )&( xi<x0(i)) );
    d1 = (x0(i) - x0(i-1));
    d2 = (xi - x0(i-1));
    x = x + (x0(i-1) + d1*f(d2/d1,alpha)).*cond;
end
cond = ( ( xi>=x0(Nlines) )&( xi<=1) );
d1 = ( 1 - x0(Nlines) );
d2 = (xi - x0(Nlines) );
x = x + ( x0(Nlines) + d1*f(d2/d1,alpha)).*cond;
% scaling
x = xmin + (xmax-xmin)*x;

end

function y = f(xi,alpha)
y  = ( exp(alpha*xi) - 1 )./( exp(alpha) - 1 );
end