function [f,df,ddf] = my_abs_fcn()
f   = @(x)   sa_f(x);
df  = @(x)  dsa_f(x);
ddf = @(x) ddsa_f(x);
end

function f = sa_f(x)
[f,~,~] = my_abs(x);
end

function df = dsa_f(x)
[~,df,~] = my_abs(x);
end

function ddf = ddsa_f(x)
[~,~,ddf] = my_abs(x);
end

function [f,df,ddf] = my_abs(x)
f  = abs(x);
df = 2*heaviside(x) - 1;
ddf = zeros(size(x));
ddf(x==0) = 2;
end