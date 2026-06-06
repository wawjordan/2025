function [f,df,ddf] = vinokur_interior_spacing_fcn(N,ti,di,refine)
si = vinokur_interior_set_spacing(N,ti,di,refine);
[ delta, xi1, b ] = vinokur_get_delta_interior(ti,si,refine);
f   = @(t)   fun(t,ti,delta,xi1,b);
df  = @(t)  dfun(t,ti,delta,xi1,b);
ddf = @(t) ddfun(t,ti,delta,xi1,b);
end
function xi = fun(t,ti,delta,xi1,b)
[xi,~,~] = vinokur_interior_f(t,ti,delta,xi1,b);
end
function dxi = dfun(t,ti,delta,xi1,b)
[~,dxi,~] = vinokur_interior_f(t,ti,delta,xi1,b);
end
function ddxi = ddfun(t,ti,delta,xi1,b)
[~,~,ddxi] = vinokur_interior_f(t,ti,delta,xi1,b);
end