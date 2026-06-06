function [f,df,ddf] = vinokur_two_sided_spacing_fcn(N,d0,d1,refine)
[s0,s1] = vinokur_two_sided_set_both_spacing(N,d0,d1,refine);
[delta,branch] = vinokur_get_delta_two_sided(s0,s1,refine);
f   = @(t)   fun(t,s0,s1,delta,branch);
df  = @(t)  dfun(t,s0,s1,delta,branch);
ddf = @(t) ddfun(t,s0,s1,delta,branch);
end
function xi = fun(t,s0,s1,delta,branch)
[xi,~,~] = vinokur_two_sided_f(t,s0,s1,delta,branch);
end
function dxi = dfun(t,s0,s1,delta,branch)
[~,dxi,~] = vinokur_two_sided_f(t,s0,s1,delta,branch);
end
function ddxi = ddfun(t,s0,s1,delta,branch)
[~,~,ddxi] = vinokur_two_sided_f(t,s0,s1,delta,branch);
end