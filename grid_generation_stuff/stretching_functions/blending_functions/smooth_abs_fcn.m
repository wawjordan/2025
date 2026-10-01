function [f,df,ddf] = smooth_abs_fcn(k)
f   = @(x)   sa_f(x,k);
df  = @(x)  dsa_f(x,k);
ddf = @(x) ddsa_f(x,k);
end

function f = sa_f(x,k)
[f,~,~] = smooth_abs(x,k);
end

function df = dsa_f(x,k)
[~,df,~] = smooth_abs(x,k);
end

function ddf = ddsa_f(x,k)
[~,~,ddf] = smooth_abs(x,k);
end