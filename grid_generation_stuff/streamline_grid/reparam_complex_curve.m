function [t,L] = reparam_complex_curve(f1,f2,N,tmin,tmax)
% creates a set of points parameterized in (approximate) arc length space
% of f1(tmin,tmax) with spacing from f2(0,1)
t0 = linspace(0,1,N);
z = f1( tmin + (tmax-tmin)*t0 ); z = z(:).';
points = [real(z);imag(z)];
tc = [ 0; cumsum( sqrt( sum( abs(points(:,2:end) - points(:,1:end-1)).^2, 1) ) ).'];
L = tc(end);
t = interp1(tc/L,t0,f2(t0));
t = tmin + (tmax-tmin)*t;
end