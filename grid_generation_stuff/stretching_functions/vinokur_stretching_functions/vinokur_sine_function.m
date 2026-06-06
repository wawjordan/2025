function [ x ] = vinokur_sine_function( y )
% Description: Approximately solves the function y = sin(x)/x for x
% if y < 0.26938972
%     x = sine_function_1(y);
% else
%     x = sine_function_2(y);
% end
mask = (y < 0.26938972);
x = zeros(size(y));
x(mask)  = sine_function_1(y(mask));
x(~mask) = sine_function_2(y(~mask));
end

function x = sine_function_1(y)
a =  1;
b = -1;
c =  1;
d = -(1 + pi^2/6);
e =  6.794732;
f = -13.205501;
g =  11.726095;
% x = pi*( a + b*y + c*y.^2 + d*y.^3 + e*y.^4 + f*y.^5 + g*y.^6 );
x = pi*(a + y.*(b + y.*(c + y.*(d + y.*(e + y.*(g*y + f))))));
end

function x = sine_function_2(y)
y1 = 1 - y;
a =  1;
b =  0.15;
c =  0.057321429;
d =  0.048774238;
e = -0.053337753;
f =  0.075845134;
% x = sqrt(6*y1).*( a + b*y1 + c*y1.^2 + d*y1.^3 + e*y1.^4 + f*y1.^5 );
x = sqrt(6*y1).*(a + y1.*(b + y1.*(c + y1.*(d + y1.*(f*y1 + e)))));
end