function x = weight_grid(xmin,xmax,N,phi)
xi = linspace(0,1,N);
% xi = linspace(-1,1,N);
xi_c = 0.5*( xi(1:N-1) + xi(2:N) );
abc = 1./phi(xi_c);
di = zeros(N-2,1);
di(1) = abc(1)*xmin;
di(N-2) = abc(N-1)*xmax;
x = zeros(N,1);
x(1) = xmin;
x(N) = xmax;
x(2:N-1) = local_tridiag(-abc(1:N-2),abc(1:N-2) + abc(2:N-1),-abc(2:N-1),di);
end

function x = local_tridiag(a,b,c,d)
N = length(b);
x = zeros(N,1);
for i = 2:N
    w = a(i)/b(i-1);
    b(i) = b(i)-w*c(i-1);
    d(i) = d(i)-w*d(i-1);
end
x(N) = d(N)/b(N);
for i = N-1:-1:1
    x(i) = (d(i)-c(i)*x(i+1))/b(i);
end
end