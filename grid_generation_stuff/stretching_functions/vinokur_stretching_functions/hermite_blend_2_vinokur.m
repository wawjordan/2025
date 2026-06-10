function [f,df,ddf] = hermite_blend_2_vinokur(N,d0,d1,off,refine)

x1 = 0.5 - off;
x2 = 0.5 + off;

% [f0,df0,ddf0] = vinokur_two_sided_spacing_fcn(N,d0,d1,refine);
if (mod(N,2)==0)
    N2 = N/2;
else
    N2 = (N-1)/2+1;
end
[f0,df0,ddf0] = vinokur_two_sided_spacing_fcn(N2,2*d0,2*d1,refine);

% hold on; fplot(@(t)0.5*f0(2*t),[0,0.5])
% hold on; fplot(@(t)1-0.5*f0(2-2*t),[0.5,1])

% f2   = @(t)1 - f0(1-t);
% df2  = @(t)df0(1-t);
% ddf2 = @(t)-ddf0(1-t);
% [f,df,ddf] = hermite_blended_functions(x1,x2,f0,f2,df0,df2,ddf0,ddf2);

f1   = @(t)0.5*f0(2*t);
df1  = @(t)df0(2*t);
ddf1 = @(t)2*ddf0(2*t);

f2   = @(t)1-0.5*f0(2-2*t);
df2  = @(t)df0(2-2*t);
ddf2 = @(t)-2*ddf0(2-2*t);
[f,df,ddf] = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end

