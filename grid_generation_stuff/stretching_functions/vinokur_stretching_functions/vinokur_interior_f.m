function [ xi, dxi, d2xi ] = vinokur_interior_f(t,ti,delta,xi1,b)
xi = xi1 + asinh( (t/ti - 1)*b )/delta;
dxi = (b/(ti*delta))./sqrt(b^2*(t/ti - 1).^2 + 1);
d2xi = -(b^3*(t/ti - 1))./(ti^2*(b^2*(t/ti - 1).^2 + 1).^(3/2))/delta;
end