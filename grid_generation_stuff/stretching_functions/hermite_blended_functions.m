function f = hermite_blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2)
f = @(t) blend(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end

function val = blend(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2)
    mask1 = t<=x1;
    mask2 = t>=x2;
    mask3 = (~mask1)&(~mask2);
    val = zeros(size(t));
    val(mask1) = f1( t(mask1) );
    val(mask2) = f2( t(mask2) );
    val(mask3) = qh_interp(t(mask3),x1,x2,f1(x1),f2(x2),df1,df2,ddf1,ddf2);
end