function x = mygeomspace( xmin, xmax, r, N )
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
gvec = r.^(0:N-2);
dx = ( xmax - xmin )/sum(gvec);
x = zeros(1,N);
x(1) = xmin;
x(N) = xmax;
for i = 2:N-1
    x(i) = x(i-1) + (r^(i-2)*dx);
end
end