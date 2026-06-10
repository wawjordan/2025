function [f,df,ddf] = hermite_blend_2_vinokur_asym(N,t0,d0,d1,off,refine)

delta_t1 = t0;
delta_t2 = 1-t0;

x1 = t0 - off*delta_t1;
x2 = t0 + off*delta_t1;

N1 = ceil(N*delta_t1);
N2 = ceil(N*delta_t2);

[f01,df01,ddf01] = vinokur_two_sided_spacing_fcn(N1,d0/delta_t1,d1/delta_t1,refine);
[f02,df02,ddf02] = vinokur_two_sided_spacing_fcn(N2,d0/delta_t2,d1/delta_t2,refine);

N_   = [N1,N2];
xi1 = [0,cumsum(N_-1)/sum(N_(:)-1)];
xi2 = [0,t0,1.0];
fh = hermite_blended_piecewise_linear(xi1,xi2,off,off);

tt1  = @(t) fh(t)/delta_t1;
f1   = @(t)       delta_t1*f01( tt1(t) );
df1  = @(t)               df01( tt1(t) );
ddf1 = @(t) (1/delta_t1)*ddf01( tt1(t) );

tt2 = @(t) ( 1 - fh(t) )/delta_t2;
f2   = @(t)     1-delta_t2*f02( tt2(t) );
df2  = @(t)               df02( tt2(t) );
ddf2 = @(t)-(1/delta_t2)*ddf02( tt2(t) );

[f,df,ddf] = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);

end
% function [f,df,ddf] = hermite_blend_2_vinokur_asym(N,t0,d0,d1,off,refine)
% 
% % x1 = 0.5 - off;
% % x2 = 0.5 + off;
% x1 = t0 - off;
% x2 = t0 + off;
% 
% delta_t1 = t0;
% delta_t2 = 1-t0;
% 
% N1 = ceil(N*delta_t1);
% % N2 = floor(N*delta_t2);
% N2 = ceil(N*delta_t2);
% 
% N_   = [N1,N2];
% xi1 = [0,cumsum(N_-1)/sum(N_(:)-1)];
% xi2 = [0,t0,1.0];
% fh = hermite_blended_piecewise_linear(xi1,xi2,off,off);
% 
% % [f0,df0,ddf0] = vinokur_two_sided_spacing_fcn(N,d0,d1,refine);
% % if (mod(N,2)==0)
% %     N2 = N/2;
% % else
% %     N2 = (N-1)/2+1;
% % end
% 
% 
% if (mod(N,2)==0)
%     N3 = N/2;
% else
%     N3 = (N-1)/2+1;
% end
% [f0,df0,ddf0] = vinokur_two_sided_spacing_fcn(N3,2*d0,2*d1,refine);
% f1_1   = @(t)0.5*f0(2*t);
% df1_1  = @(t)df0(2*t);
% ddf1_1 = @(t)2*ddf0(2*t);
% 
% f2_1   = @(t)1-0.5*f0(2-2*t);
% df2_1  = @(t)df0(2-2*t);
% ddf2_1 = @(t)-2*ddf0(2-2*t);
% 
% [f01,df01,ddf01] = vinokur_two_sided_spacing_fcn(N1,d0/delta_t1,d1/delta_t1,refine);
% 
% 
% tt1 = @(t) fh(t)/delta_t1;
% f1   = @(t) delta_t1*f01( tt1(t) );
% df1  = @(t) df01( tt1(t) );
% ddf1 = @(t) (1/delta_t1)*ddf01( tt1(t) );
% 
% 
% [f02,df02,ddf02] = vinokur_two_sided_spacing_fcn(N2,d0/delta_t2,d1/delta_t2,refine);
% 
% tt2 = @(t) ( 1 - fh(t) )/delta_t2;
% f2   = @(t) 1-delta_t2*f02( tt2(t) );
% df2  = @(t) df02( tt2(t) );
% ddf2 = @(t)-(1/delta_t2)*ddf02( tt2(t) );
% [f,df,ddf] = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2);
% end

