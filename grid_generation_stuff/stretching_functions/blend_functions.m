function f = blend_functions(x1,x2,f1,f2)
f = @(t) blend(t,x1,x2,f1,f2);
end

function val = blend(t,x1,x2,f1,f2)
    mask1 = t<=x1;
    mask2 = t>=x2;
    mask3 = (~mask1)&(~mask2);
    val = zeros(size(t));
    val(mask1) = f1( t(mask1) );
    val(mask2) = f2( t(mask2) );
    tt = ( t(mask3) - x1 )/(x2-x1);
    pt = tt.^3.*(tt.*(  6*tt -  15) + 10);
    val(mask3) = (1-pt) .* f1( t(mask3) ) + pt .* f2( t(mask3) );
end