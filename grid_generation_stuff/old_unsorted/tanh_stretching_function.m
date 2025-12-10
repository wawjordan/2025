function F = tanh_stretching_function(x1,x2,d1,d2,N)
% x1 and x2 are the start/end points
% d1 and d2 are spacings at the endpoints
% N is number of points

L = abs(x2-x1);


xmap = @(x) (x-x1)/L;

xmapinv = @(x) L*x + x1;

% get approximate (normalized) slopes
% S0 = d1*(N-1)/L;
% S1 = d2*(N-1)/L;

S0 = L/(d1*(N-1));
S1 = L/(d2*(N-1));

% get A and B
A = sqrt(S0/S1);
B = sqrt(S0*S1);

small = 1e-12;

if B > 1 + small
    dy = invB1(B);
    u = @(x) 0.5*( 1 + tanh(dy*(x-0.5))/tanh(0.5*dy) );
elseif B < 1 - small
    dy = invB2(B);
    u = @(x) 0.5*( 1 + tan(dy*(x-0.5))/tan(0.5*dy) );
else
    u = @(x) x.*(1 + 2*(B-1).*(x-0.5).*(1-x));
end

t = @(x) u(x)./( A + (1-A)*u(x) );

F = @(x) xmapinv(t(xmap(x)));

end


function x = invB1(y)

% get guess
a = 2.7829681;
if y < a
    x0 = invB1_a(y);
else
    x0 = invB1_b(y);
end
x = x0;
% x = fzero(@(z)sinh(z)./z - y,x0);

end

function x = invB1_a(y)
y1 = y - 1;
a =  1;
b = -0.15;
c =  0.057321429;
d = -0.024907295;
e =  0.0077424461;
f = -0.0010794123;

x = sqrt(6*y1).*( a + b*y1 + c*y1.^2 + d*y1.^3 + e*y1.^4 + f*y1.^5 );
end

function x = invB1_b(y)
v = log(y);
w = 1./y - 0.028527431;
a = -0.02041793;
b =  0.24902722;
c =  1.9496443;
d = -2.6294547;
e =  8.56795911;

x = v + ( 1 + 1./v ).*log(2*v) + a + b*w + c*w.^2 + d*w.^3 + e*w.^4;
end

function x = invB2(y)
a = 0.26938972;
if y < a
    x0 = invB2_a(y);
else
    x0 = invB2_b(y);
end
x = x0;
% x = fzero(@(z)sin(z)./z - y,x0);
end

function x = invB2_a(y)
a =  1;
b = -1;
c =  1;
d = -(1 + pi^2/6);
e =  6.794732;
f = -13.205501;
g =  11.726095;

x = pi*( a + b*y + c*y.^2 + d*y.^3 + e*y.^4 + f*y.^5 + g*y.^6 );
end


function x = invB2_b(y)
y1 = y - 1;
a =  1;
b =  0.15;
c =  0.057321429;
d =  0.048774238;
e = -0.053337753;
f =  0.075845134;

x = sqrt(6*y1).*( a + b*y1 + c*y1.^2 + d*y1.^3 + e*y1.^4 + f*y1.^5 );
end