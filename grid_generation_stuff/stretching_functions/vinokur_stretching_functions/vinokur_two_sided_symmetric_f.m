function [xi,dxi,ddxi,s] = vinokur_two_sided_symmetric_f(t,s,N)
if nargin==3
    delta_xi = s;
    s = vinokur_two_sided_symmetric_set_mid_spacing(N,delta_xi);
end
[xi,dxi,ddxi] = vinokur_two_sided_symmetric(t,s);
end

function [s,d] = vinokur_two_sided_symmetric_set_mid_spacing(N,delta_xi)
if ~( mod(N,1)<eps(1) && N>1 ), error('N must be a postive integer >2'); end
dt = 1/(N-1);    
if N<4
    warning('N must be >4 to determine s; returning s=1 (linear)');
    s = 1;
    if nargout > 1
        d = dt;
    end
    return
end
if mod(N,2)==0 % even number of edges --> odd number of bins
    if (delta_xi>1)
        error('delta_xi too big')
    end
    t0 = 0.5-dt/2;
    t1 = 0.5+dt/2;
    % initial guesses for slope
    % space left for N/2-1 points on 1 side = (1-dx)/2
    % 
    s00 = ( N/2 - 1 ) / ( (1-delta_xi)/2 );
else
    if (delta_xi>0.5)
        error('delta_xi too big')
    end
    t0 = 0.5-dt;
    t1 = 0.5;
    % initial guesses for slope
    % space left for (N-1)/2 -1 points on 1 side = (1-2*dx)/2
    % 
    s00 = ( (N-1)/2 - 1 ) / ( (1-2*delta_xi)/2 );
end
s0 = s00/10000;
s1 = s00*10000;


options = optimset();
s = fzero( @(s) obj_fun(s,t0,t1,delta_xi),[s0,s1], options );
if nargout > 1
    [~,d] = obj_fun(s,t0,t1,delta_xi);
end
    function [e,d] = obj_fun(s,t0,t1,delta_target)
        xi = vinokur_two_sided_symmetric([t0,t1],s);
        d = xi(2)-xi(1);
        e = d - delta_target;
    end
end

function [xi,dxi,ddxi] = vinokur_two_sided_symmetric(t,s)
tol = 0.001;
B = 1/s;
if abs(B-1) < tol
    xi = t.*(1 + 2*(B-1)*(t-0.5).*(1-t));
    if nargout>1
        dxi  = 6*(B-1)*t.*(1-t) - B + 2;
        if nargout > 2
            ddxi = 6*(B-1)*(1-2*t);
        end
    end
elseif B > 1
    delta = sinh_function( B );
    xi   = 0.5*( 1 + tanh(delta*(t-0.5)) / tanh(0.5*delta) );
    if nargout>1
        dxi  = 0.5*((delta*( 1 - tanh(delta*(t-0.5) ).^2 ) ) ...
             / tanh(0.5*delta));
        if nargout > 2
            ddxi = ( delta^2*tanh(delta*(t-0.5)).* ...
                     (tanh(delta*(t-0.5)).^2-1) ) / tanh(0.5*delta);
        end
    end
elseif B < 1
    delta = sine_function( B );
    xi = 0.5*( 1 + tan(delta*(t-0.5)) / tan(0.5*delta) );
    if nargout>1
        dxi = 0.5*( (delta*(1 + tan(delta*(t-0.5)).^2) ) / tan(0.5*delta) );
        if nargout > 2
            ddxi = ( delta^2*tan(delta*(t-0.5)).* ...
                     (tan(delta*(t-0.5)).^2+1)) / tan(0.5*delta);
        end
    end
end
end


function x = sine_function( y )
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

function x = sinh_function( y )
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