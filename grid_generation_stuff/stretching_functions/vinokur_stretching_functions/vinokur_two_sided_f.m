function [xi,dxi,d2xi] = vinokur_two_sided_f(t,s0,s1,delta,branch)
s0_ = 1/s0;
s1_ = 1/s1;
A   = sqrt(s0_/s1_);
B   = sqrt(s0_*s1_);
switch branch
    case 0
        u = t.*(1 + 2*(B-1)*(t-0.5).*(1-t));
        du = 6*(B-1)*t.*(1-t) - B + 2;
        du2 = 6*(B-1)*(1-2*t);
    case 1
        u  = 0.5*( 1 + tanh(delta*(t-0.5))/tanh(0.5*delta) );
        du = 0.5*((delta*( 1 - tanh(delta*(t-0.5) ).^2 ) )/tanh(0.5*delta));
        du2 = (delta^2*tanh(delta*(t-0.5)).*(tanh(delta*(t-0.5)).^2-1))/tanh(0.5*delta);
    case 2
        u  = 0.5*( 1 + tan(delta*(t-0.5))/tan(0.5*delta) );
        du = 0.5*((delta*(1 + tan(delta*(t-0.5)).^2))/tan(0.5*delta));
        du2 = (delta^2*tan(delta*(t-0.5)).*(tan(delta*(t-0.5)).^2+1))/tan(0.5*delta);
end
xi  = u./( A + (1-A)*u );
dxi = du./(A - u*(A - 1)) + ((A-1)*u.*du)./(A - (A-1)*u).^2;
d2 = (2*(A-1))./(A - u*(A - 1)).^2 + (2*u*(A - 1)^2)./(A - u*(A - 1)).^3;
d2xi = du.^2 .* d2 + du2.*dxi;
end