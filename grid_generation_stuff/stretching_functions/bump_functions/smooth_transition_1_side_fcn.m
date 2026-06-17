function [f,df,ddf] = smooth_transition_1_side_fcn(a,b,h1,h2)
f   = @(t)   fun(t,a,b,h1,h2);
df  = @(t)  dfun(t,a,b,h1,h2);
ddf = @(t) ddfun(t,a,b,h1,h2);
end
function xi = fun(t,a,b,h1,h2)
[xi,~,~] = smooth_transition_1_side(t,a,b,h1,h2);
end
function dxi = dfun(t,a,b,h1,h2)
[~,dxi,~] = smooth_transition_1_side(t,a,b,h1,h2);
end
function ddxi = ddfun(t,a,b,h1,h2)
[~,~,ddxi] = smooth_transition_1_side(t,a,b,h1,h2);
end