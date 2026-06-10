function [f,df,ddf] = hermite_blend_2_vinokur_asym(N,t0,d0,d1,off,refine)

x1 = t0 - off * 2*t0;
x2 = t0 + off * 2*(1-t0);

N1 = ceil(N*t0);
N2 = floor(N*(1-t0));

% N = [N1,N2];
% x = [0,cumsum(N-1)/sum(N(:)-1)];
% y = [0,t0,1];
% dm = off;
% dp = off;
% [f_l,df_l,ddf_l] = hermite_blended_piecewise_linear(x,y,dm,dp);

[f01,df01,ddf01] = vinokur_two_sided_spacing_fcn(N1,2*d0,2*d1,refine);
[f02,df02,ddf02] = vinokur_two_sided_spacing_fcn(N2,2*d0,2*d1,refine);

% f1   = @(t)0.5*f01(2*t);
% df1  = @(t)df01(2*t);
% ddf1 = @(t)2*ddf01(2*t);
% 
% f2   = @(t)1-0.5*f02(2-2*t);
% df2  = @(t)df02(2-2*t);
% ddf2 = @(t)-2*ddf02(2-2*t);


f1   = @(t)t0*f01(2*t);
df1  = @(t)2*t0*df01(2*t);
ddf1 = @(t)4*t0*ddf01(2*t);

f2   = @(t)1-(1-t0)*f02(2-2*t);
df2  = @(t)2*(1-t0)*df02(2-2*t);
ddf2 = @(t)-4*(1-t0)*ddf02(2-2*t);

[f,df,ddf] = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end

