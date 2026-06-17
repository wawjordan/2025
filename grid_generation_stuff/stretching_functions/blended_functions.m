function [f,df,ddf] = blended_functions(x1,x2,f1,f2,df1,df2,ddf1,ddf2)
f   = @(t) fun(  t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
df  = @(t) dfun( t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
ddf = @(t) ddfun(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end
function fval = fun(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2)
[fval,~,~] = blend(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end
function dfval = dfun(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2)
[~,dfval,~] = blend(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end
function ddfval = ddfun(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2)
[~,~,ddfval] = blend(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2);
end
function [val,dval,ddval] = blend(t,x1,x2,f1,f2,df1,df2,ddf1,ddf2)
    mask1 = t<=x1;
    mask2 = t>=x2;
    mask3 = (~mask1)&(~mask2);
    val   = zeros(size(t));
    dval  = zeros(size(t));
    ddval = zeros(size(t));
    val(mask1)   = f1( t(mask1) );
    dval(mask1)  = df1( t(mask1) );
    ddval(mask1) = ddf1( t(mask1) );
    val(mask2)   = f2( t(mask2) );
    dval(mask2)  = df2( t(mask2) );
    ddval(mask2) = ddf2( t(mask2) );
    % [val(mask3),dval(mask3),ddval(mask3)] = qh_interp(t(mask3),x1,x2,f1(x1),f2(x2),df1(x1),df2(x2),ddf1(x1),ddf2(x2));
    [fb,dfb,ddfb] = smooth_transition_1_side(t(mask3),x1,x2,1,0);
    val(mask3) = f1( t(mask3) ).*        fb ...
               + f2( t(mask3) ).* ( 1 -  fb );
    dval(mask3) = df1( t(mask3) ).*      fb ...
                + df2( t(mask3) ).*( 1 - fb ) ...
                + ( f1( t(mask3) ) -  f2( t(mask3) ) ).* dfb;
   ddval(mask3) = ddf1( t(mask3) ).*      fb ...
                +2*df1( t(mask3) ).*     dfb ...
                +   f1( t(mask3) ).*    ddfb ...
                + ddf2( t(mask3) ).*( 1 - fb ) ...
                -2*df2( t(mask3) ).*     dfb ...
                -   f2( t(mask3) ).*    ddfb;
end