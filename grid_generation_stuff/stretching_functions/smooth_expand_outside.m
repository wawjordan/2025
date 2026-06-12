function [f,df] = smooth_expand_outside(t,N1,N2,a,b,c,d,dh1,dh2,h0,h1,h2,refine)
% fplot(@(t)f(2*(t-0.5)),[0.5,1])
% hold on
% fplot(@(t)f(1-2*t),[0,0.5])

s1               = vinokur_one_sided_set_spacing(N1,dh1/(h1-h0),refine);
[delta1,branch1] = vinokur_get_delta_one_sided(s1,refine);
s2               = vinokur_one_sided_set_spacing(N2,dh2/(h2-h0),refine);
[delta2,branch2] = vinokur_get_delta_one_sided(s2,refine);
f = zeros(size(t));
m1 = (t<b);
m2 = (t>=b)&(t<=c);
m3 = (t>c);
% [f1,~,~] = vinokur_one_sided_f(1-2*t(m1),s1,delta1,branch1);
% [f2,~,~] = smooth_transition_1_side(t(m1),a,b,0,1);
% [f3,~,~] = smooth_transition_1_side(t(m3),c,d,1,0);
% [f4,~,~] = vinokur_one_sided_f(2*t(m3)-1,s2,delta2,branch2);
% f(m1) = ( h0 + (h1-h0)*f1 ) .* ( 1 - f2 ) + h0*f2;
% f(m2) = h0;
% f(m3) = ( h0 + (h2-h0)*f4 ) .* ( 1 - f3 ) + h0*f3;

[f1,df1,ddf1] = vinokur_one_sided_f((b-t(m1))/b,s1,delta1,branch1);
[f2,df2,ddf2] = smooth_transition_1_side(t(m1),a,b,0,1);
[f3,df3,ddf3] = smooth_transition_1_side(t(m3),c,d,1,0);
[f4,df4,ddf4] = vinokur_one_sided_f((t(m3)-c)/(1-c),s2,delta2,branch2);
% f1 = (a-t(m1))/a;
% f4 = (t(m3)-d)/(1-d);
f(m1) = ( h0 + (h1-h0)*f1 ) .* ( 1 - f2 ) + h0*f2;
f(m2) = h0;
f(m3) = ( h0 + (h2-h0)*f4 ) .* ( 1 - f3 ) + h0*f3;

df(m1) = (h1-h0)*(-1/a)*df1 .* ( 1 - f2 ) - ( h0 + (h1-h0)*f1 ) .* df2 + h0*df2;
df(m2) = 0;
df(m3) = (h2-h0)*(1/(1-d))*df4 .* (1 - f3 ) - ( h0 + (h2-h0)*f4 ) .* df3 + h0*df3;
% f(m1) = f2;
% f(m2) = h0;
% f(m3) = f3;
end