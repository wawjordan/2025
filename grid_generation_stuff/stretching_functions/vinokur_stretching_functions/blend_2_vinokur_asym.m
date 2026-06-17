function [f,df,ddf] = blend_2_vinokur_asym(N,t0,d0,d1,off,refine)

delta_t1 = t0;
delta_t2 = 1-t0;

x10 = max(t0  - off*delta_t1,0);
x1  = max(x10 - 9*off*delta_t1,0);
x20 = min(t0  + off*delta_t2,1);
x2  = min(x20 + 9*off*delta_t2,1);



N1 = ceil(N*delta_t1);
N2 = ceil(N*delta_t2);

[f01,df01,ddf01] = vinokur_two_sided_spacing_fcn(N1,d0/delta_t1,d1/delta_t1,refine);
[f02,df02,ddf02] = vinokur_two_sided_spacing_fcn(N2,d0/delta_t2,d1/delta_t2,refine);

N_   = [N1,N2];
xi1 = [0,cumsum(N_-1)/sum(N_(:)-1)];
xi2 = [0,t0,1.0];
[fh,dfh,ddfh] = hermite_blended_piecewise_linear(xi1,xi2,off,off);

tt1   = @(t)   fh(t)/delta_t1;
dtt1  = @(t)  dfh(t)/delta_t1;
ddtt1 = @(t) ddfh(t)/delta_t1;
f1   = @(t)            delta_t1*f01( tt1(t) );
df1  = @(t)            delta_t1*df01( tt1(t) ).*dtt1(t);
ddf1 = @(t) delta_t1*( dtt1(t).^2 .* ddf01( tt1(t) ) + ddtt1(t).*df01(tt1(t)) );

tt2   = @(t) ( 1 - fh(t) )/delta_t2;
dtt2  = @(t)     -  dfh(t)/delta_t2;
ddtt2 = @(t)     - ddfh(t)/delta_t2;
f2   = @(t)     1-delta_t2*f02( tt2(t) );
df2  = @(t)               -delta_t2*( df02( tt2(t) ).*dtt2(t) );
ddf2 = @(t) - delta_t2 * ( dtt2(t).^2 .* ddf02( tt2(t) ) + ddtt2(t).*df02(tt2(t)) );

fmid   = @(t) f1(t0) + ones(size(t))*0.5*(f1(t0)+f2(t0));
dfmid  = @(t) zeros(size(t));
ddfmid = @(t) zeros(size(t));
[fL,dfL,ddfL] = blended_functions(x1,x10,f1,fmid,df1,dfmid,ddf1,ddfmid);
[f,df,ddf]    = blended_functions(x20,x2,fL,f2,dfL,df2,ddfL,ddf2);

% [f,df,ddf] = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);

end