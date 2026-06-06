function [f,df,ddf] = c2_blend_poly_f()
f   = @(t)   fun(t);
df  = @(t)  dfun(t);
ddf = @(t) ddfun(t);
end

function f = fun(t)
    [f,~,~] = c2_blend_poly(t);
end
function df = dfun(t)
    [~,df,~] = c2_blend_poly(t);
end
function ddf = ddfun(t)
    [~,~,ddf] = c2_blend_poly(t);
end