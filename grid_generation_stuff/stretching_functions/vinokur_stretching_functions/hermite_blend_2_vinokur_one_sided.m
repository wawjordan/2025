function [f,df,ddf] = hermite_blend_2_vinokur_one_sided(N,d0,off,refine)

d0 = get_optimized_spacing(N,d0,off,refine);

x1 = 0.5 - off;
x2 = 0.5 + off;
[f1,f2,df1,df2,ddf1,ddf2] = get_half_functions(N,d0,refine);
[f,df,ddf] = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);
% [f,df,ddf] = blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end

function [f1,f2,df1,df2,ddf1,ddf2] = get_half_functions(N,d0,refine)
if (mod(N,2)==0)
    N2 = N/2;
else
    N2 = (N-1)/2+1;
end
[f0,df0,ddf0] = vinokur_one_sided_spacing_fcn(N2,2*d0,refine);
f1   = @(t)0.5-0.5*f0(1-2*t);
df1  = @(t)       df0(1-2*t);
ddf1 = @(t)   -2*ddf0(1-2*t);
f2   = @(t)0.5+0.5*f0(2*t-1);
df2  = @(t)       df0(2*t-1);
ddf2 = @(t)    2*ddf0(2*t-1);
end

function d0 = get_optimized_spacing(N,d0_target,off,refine)
d0 = fzero( @(d)eval_spacing(N,d,off,refine)-d0_target,d0_target);
end

function d = eval_spacing(N,d,off,refine)
[f1,f2,df1,df2,ddf1,ddf2] = get_half_functions(N,d,refine);
d = blend_deriv(0.5,0.5-off,0.5+off,f1,f2,df1,df2,ddf1,ddf2);
d = d/(N-1);
end

function dval = blend_deriv(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2)
    mask1 = t<=x1;
    mask2 = t>=x2;
    mask3 = (~mask1)&(~mask2);
    dval  = zeros(size(t));
    dval(mask1)  = df1( t(mask1) );
    dval(mask2)  = df2( t(mask2) );
    [~,dval(mask3),~] = qh_interp(t(mask3),x1,x2,f1(x1),f2(x2),df1(x1),df2(x2),ddf1(x1),ddf2(x2));
end

% function d = eval_spacing(N,d,off,refine)
% [f1,f2,df1,df2,ddf1,ddf2] = get_half_functions(N,d,refine);
% [~,df,~] = blended_functions(0.5-off,0.5+off,f1,f2,df1,df2,ddf1,ddf2);
% d = df(0.5)/(N-1);
% end
