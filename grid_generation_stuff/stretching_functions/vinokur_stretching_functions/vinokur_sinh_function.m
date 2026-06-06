function [ x ] = vinokur_sinh_function( y )
% Description: Approximately solves the function y = sinh(x)/x for x
% if y < 2.7829681
%     x = sinh_function_1(y);
% else
%     x = sinh_function_2(y);
% end
mask = (y < 2.7829681);
x = zeros(size(y));
x(mask)  = sinh_function_1(y(mask));
x(~mask) = sinh_function_2(y(~mask));
end

function x = sinh_function_1(y)
y1 = y - 1;
a =  1;
b = -0.15;
c =  0.057321429;
d = -0.024907295;
e =  0.0077424461;
f = -0.0010794123;
% x = sqrt(6*y1).*( a + b*y1 + c*y1.^2 + d*y1.^3 + e*y1.^4 + f*y1.^5 );
x  = sqrt(6*y1).*(a + y1*(b + y1*(c + y1*(d + y1*(f*y1 + e)))));
end

function x = sinh_function_2(y)
v = log(y);
w = 1./y - 0.028527431;
a = -0.02041793;
b =  0.24902722;
c =  1.9496443;
d = -2.6294547;
e =  8.56795911;
% x = v + ( 1 + 1./v ).*log(2*v) + a + b*w + c*w.^2 + d*w.^3 + e*w.^4;
x = v + ( 1 + 1./v ).*log(2*v) + a + w.*(b + w.*(c + w.*(e*w + d)));
end