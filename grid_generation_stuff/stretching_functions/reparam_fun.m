function [t,L] = reparam_fun(dim,f1,f2,N,tmin,tmax)
% creates a set of points parameterized in arc length space of f1 with spacing from f2
t0 = linspace(tmin,tmax,N);
points = zeros(dim,N);
switch dim
    case(1)
        points = f1(t0);
    case(2)
        [points(1,:),points(2,:)] = f1(t0);
    case(3)
        [points(1,:),points(2,:),points(2,:)] = f1(t0);
end
% t1 = centri_param(points,mu);
t1 = [ 0; cumsum( sqrt( sum( abs(points(:,2:end) - points(:,1:end-1)).^2, 1) ) ).'];
L = t1(end);
t = interp1(t1/L,linspace(0,1,N),f2(linspace(0,1,N)));
end
