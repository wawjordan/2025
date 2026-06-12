function [f,df,ddf] = vinokur_one_sided_spacing_fcn(N,d0,refine)
s0 = vinokur_one_sided_set_spacing(N,d0,refine);
[delta,branch] = vinokur_get_delta_one_sided(s0,refine);
f   = @(t)   fun(t,s0,delta,branch);
df  = @(t)  dfun(t,s0,delta,branch);
ddf = @(t) ddfun(t,s0,delta,branch);
end
function xi = fun(t,s0,delta,branch)
[xi,~,~] = vinokur_one_sided_f(t,s0,delta,branch);
end
function dxi = dfun(t,s0,delta,branch)
[~,dxi,~] = vinokur_one_sided_f(t,s0,delta,branch);
end
function ddxi = ddfun(t,s0,delta,branch)
[~,~,ddxi] = vinokur_one_sided_f(t,s0,delta,branch);
end