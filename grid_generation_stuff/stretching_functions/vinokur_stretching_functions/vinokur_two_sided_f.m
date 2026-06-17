function [xi,dxi,d2xi] = vinokur_two_sided_f(t,s0,s1,delta,branch)
s0_ = 1/s0;
s1_ = 1/s1;
A   = sqrt(s0_/s1_);
B   = sqrt(s0_*s1_);
switch branch
    case 0
        u       = t.*(1 + 2*(B-1)*(t-0.5).*(1-t));
        du_dt   = 6*(B-1)*t.*(1-t) - B + 2;
        du2_dt2 = 6*(B-1)*(1-2*t);
    case 1
        u       = 0.5*( 1 + tanh(delta*(t-0.5))/tanh(0.5*delta) );
        du_dt   = 0.5*((delta*( 1 - tanh(delta*(t-0.5) ).^2 ) )/tanh(0.5*delta));
        du2_dt2 = (delta^2*tanh(delta*(t-0.5)).*(tanh(delta*(t-0.5)).^2-1))/tanh(0.5*delta);
    case 2
        u       = 0.5*( 1 + tan(delta*(t-0.5))/tan(0.5*delta) );
        du_dt   = 0.5*((delta*(1 + tan(delta*(t-0.5)).^2))/tan(0.5*delta));
        du2_dt2 = (delta^2*tan(delta*(t-0.5)).*(tan(delta*(t-0.5)).^2+1))/tan(0.5*delta);
end
xi    = u./( A + (1-A)*u );
dxi_du = 1./(A - u*(A - 1)) + (u*(A - 1))./(A - u*(A - 1)).^2;
% dxi = du./(A - u*(A - 1)) + ((A-1)*u.*du)./(A - (A-1)*u).^2;
dxi = dxi_du.*du_dt;
d2xi_du2 = (2*(A-1))./(A - u*(A - 1)).^2 + (2*u*(A - 1)^2)./(A - u*(A - 1)).^3;
d2xi = du_dt.^2 .* d2xi_du2 + du2_dt2.*dxi_du;
end