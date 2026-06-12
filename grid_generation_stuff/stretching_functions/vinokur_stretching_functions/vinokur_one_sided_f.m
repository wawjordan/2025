function [xi,dxi,d2xi] = vinokur_one_sided_f(t,s0,delta,branch)
s0_ = 1/s0;
switch branch
    case 0
        xi  = t.*(1-0.5*(s0_-1)*(1-t).*(2-t));
        dxi = 1 - (t*(1.5*t - 3) + 1)*(s0_-1);
        d2xi = -(3*t-3)*(s0_-1);
    case 1
        xi = 1 + ( tanh(delta*(t-1))/tanh(delta) );
        dxi = (delta/tanh(delta))*(1 - tanh(delta*(t-1)).^2);
        d2xi = (2*delta^2*tanh(delta*(t-1)).*(tanh(delta*(t-1)).^2 - 1))/tanh(delta);
    case 2
        xi = 1 + ( tan(delta*(t-1))/tan(delta) );
        dxi = (delta/tan(delta))*(1 + tan(delta*(t-1)).^2);
        d2xi = (2*delta^2*tan(delta*(t-1)).*(tan(delta*(t-1)).^2 + 1))/tan(delta);
end
end