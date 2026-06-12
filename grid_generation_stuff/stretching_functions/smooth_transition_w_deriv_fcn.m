function [f,df,ddf] = smooth_transition_w_deriv_fcn(h1,h2,a,b,c,d)
f   = @(t)   fun(t,h1,h2,a,b,c,d);
df  = @(t)  dfun(t,h1,h2,a,b,c,d);
ddf = @(t) ddfun(t,h1,h2,a,b,c,d);
end
function val = fun(t,h1,h2,a,b,c,d)
[val,~,~] = smooth_transition_w_deriv(t,h1,h2,a,b,c,d);
end
function val = dfun(t,h1,h2,a,b,c,d)
[~,val,~] = smooth_transition_w_deriv(t,h1,h2,a,b,c,d);
end
function val = ddfun(t,h1,h2,a,b,c,d)
[~,~,val] = smooth_transition_w_deriv(t,h1,h2,a,b,c,d);
end