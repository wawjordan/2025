function [f,df,ddf] = smoothmax_smu_fcn(f1,f2,df1,df2,ddf1,ddf2,e)
f   = @(x)   sm_f(x,f1,f2,df1,df2,ddf1,ddf2,e);
df  = @(x)  dsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e);
ddf = @(x) ddsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e);
end

function f = sm_f(x,f1,f2,df1,df2,ddf1,ddf2,e)
[f,~,~] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);
end

function df = dsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e)
[~,df,~] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);
end

function ddf = ddsm_f(x,f1,f2,df1,df2,ddf1,ddf2,e)
[~,~,ddf] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e);
end