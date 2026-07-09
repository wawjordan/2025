function [f,df,ddf] = smooth_ramps(t,xb,vals)
f   = zeros(size(t));
df  = zeros(size(t));
ddf = zeros(size(t));
n_vals = numel(vals);
mask = t<xb(1);
f(mask) = vals(1);
for i = 2:n_vals
    mask = (t>=xb(i-1))&(t<=xb(i));
    [f1,df1,ddf1] = smooth_transition_1_side(t(mask),xb(i-1),xb(i),vals(i-1),vals(i));
    f(mask) = f1;
    df(mask) = df1;
    ddf(mask) = ddf1;
end
mask = t>xb(n_vals);
f(mask) = vals(n_vals);
end