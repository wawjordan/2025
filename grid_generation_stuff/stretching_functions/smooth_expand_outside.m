function [f,df] = smooth_expand_outside(t,N1,N2,a,b,c,d,dh1,dh2,h0,h1,h2,refine)

% s1               = vinokur_one_sided_set_spacing(N1,dh1/(h1-h0),refine);
s1 = 1000000;
[delta1,branch1] = vinokur_get_delta_one_sided(s1,refine);
% s2               = vinokur_one_sided_set_spacing(N2,dh2/(h2-h0),refine);
s2 = 1000000;
[delta2,branch2] = vinokur_get_delta_one_sided(s2,refine);

m1 = (t<b);
m2 = (t>=b)&(t<=c);
m3 = (t>c);
[f1,df1,~] = vinokur_one_sided_f((b-t(m1))/b,s1,delta1,branch1);
[f2,df2,~] = smooth_transition_1_side(t(m1),a,b,0,1);
[f3,df3,~] = smooth_transition_1_side(t(m3),c,d,1,0);
[f4,df4,~] = vinokur_one_sided_f((t(m3)-c)/(1-c),s2,delta2,branch2);

f = zeros(size(t));
f(m1) = ( h0 + (h1-h0)*f1 ) .* ( 1 - f2 ) + h0*f2;
f(m2) = h0;
f(m3) = ( h0 + (h2-h0)*f4 ) .* ( 1 - f3 ) + h0*f3;

df(m1) = (h1-h0)*(-1/a)*df1 .* ( 1 - f2 ) - ( h0 + (h1-h0)*f1 ) .* df2 + h0*df2;
df(m2) = 0;
df(m3) = (h2-h0)*(1/(1-d))*df4 .* (1 - f3 ) - ( h0 + (h2-h0)*f4 ) .* df3 + h0*df3;
end