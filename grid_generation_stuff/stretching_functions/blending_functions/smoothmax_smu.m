function [f,df,ddf] = smoothmax_smu(x,f1,f2,df1,df2,ddf1,ddf2,e)
f1x = f1(x);
f2x = f2(x);
df1x = df1(x);
df2x = df2(x);

ddf1x = ddf1(x);
ddf2x = ddf2(x);

f   = 0.5*( f1x + f2x + sqrt( (f1x-f2x).^2 + e ) );
df  = 0.5*(df1x + df2x + ((df1x - df2x ).*(f1x - f2x))./sqrt( (f1x - f2x).^2 + e ) );
ddf = 0.5*( df1x - df2x ).^2 ./ sqrt( (f1x - f2x).^2 + e ) ...
    + 0.5*( ddf1x + ddf2x ) ...
    - 0.5*( (df1x - df2x).^2 .* (f1x - f2x).^2 )./(sqrt( (f1x - f2x).^2 + e )).^3 ...
    + 0.5*( (ddf1x  - ddf2x) .* (f1x - f2x) )./sqrt( (f1x - f2x).^2 + e );

end
 % f = ( a + b + sqrt( (a-b)^2 + e ) )/2
 % diff(a(x), x)/2 + diff(b(x), x)/2 + ((diff(a(x), x) - diff(b(x), x))*(a(x) - b(x)))/(2*(e + (a(x) - b(x))^2)^(1/2))
 
 %  (diff(a(x), x) - diff(b(x), x))^2/(2*(e + (a(x) - b(x))^2)^(1/2)) ...
 % + diff(a(x), x, x)/2 ...
 % + diff(b(x), x, x)/2 ...
 % - ((diff(a(x), x) - diff(b(x), x))^2*(a(x) - b(x))^2)/(2*(e + (a(x) - b(x))^2)^(3/2)) ...
 % + ((diff(a(x), x, x) - diff(b(x), x, x))*(a(x) - b(x)))/(2*(e + (a(x) - b(x))^2)^(1/2))