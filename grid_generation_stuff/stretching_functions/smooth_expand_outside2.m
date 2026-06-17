function [f,df] = smooth_expand_outside2(t,N1,N2,a,b,c,d,dh1,dh2,h0,h1,h2)

% s1               = vinokur_one_sided_set_spacing(N1,dh1/(h1-h0),refine);
% [delta1,branch1] = vinokur_get_delta_one_sided(s1,refine);
% s2               = vinokur_one_sided_set_spacing(N2,dh2/(h2-h0),refine);
% [delta2,branch2] = vinokur_get_delta_one_sided(s2,refine);



m1 = (t<b);
m2 = (t>=b)&(t<=c);
m3 = (t>c);
% [f1,df1,~] = vinokur_one_sided_f((b-t(m1))/b,s1,delta1,branch1);
% [~,~,dx2] = my_geomspace(N1,0,dx0=dh1,xmax=h1);
% [s0,s1] = vinokur_two_sided_set_both_spacing(N1,dh1,dx2,refine);
% [delta,branch] = vinokur_get_delta_two_sided(s0,s1,refine);
% [f1,df1,~] = vinokur_two_sided_f((b-t(m1))/b,s0,s1,delta,branch);

% [~,~,~,f01] = my_geomspace_w_refine(N1,1,0,xmax=h1,dx0=dh1);
% rN = r^(N-1);
% val = xmin + dx0*(rN.^x - 1)/(r-1);
[~,r1,~] = my_geomspace(N1,h0,xmax=h1,dx0=dh1);
f1 = h0 + dh1*( (r1^(N1-1) ).^( (b-t(m1))/b ) )/(r1-1);

[~,r4,~] = my_geomspace(N2,h0,xmax=h2,dx0=dh2);
f4 = h0 + dh2*( (r4^(N2-1) ).^( (t(m3)-c)/(1-c) ) )/(r4-1);

[f2,df2,~] = smooth_transition_1_side(t(m1),a,b,0,1);
[f3,df3,~] = smooth_transition_1_side(t(m3),c,d,1,0);


f = zeros(size(t));
df = zeros(size(t));
f(m1) = ( h0 + (h1-h0)*f1 ) .* ( 1 - f2 ) + h0*f2;
f(m2) = h0;
f(m3) = ( h0 + (h2-h0)*f4 ) .* ( 1 - f3 ) + h0*f3;

% df(m1) = (h1-h0)*(-1/a)*df1 .* ( 1 - f2 ) - ( h0 + (h1-h0)*f1 ) .* df2 + h0*df2;
% df(m2) = 0;
% df(m3) = (h2-h0)*(1/(1-d))*df4 .* (1 - f3 ) - ( h0 + (h2-h0)*f4 ) .* df3 + h0*df3;
end