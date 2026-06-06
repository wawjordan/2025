function val = hermite_blend_segments(t,x1,x2,y1,y2,s1,s2)
    mask1 = t<=x1;
    mask2 = t>=x2;
    mask3 = (~mask1)&(~mask2);
    val = zeros(size(t));
    val(mask1) = y1 + s1*(t(mask1)-x1);
    val(mask2) = y2 + s2*(t(mask2)-x2);
    val(mask3) = qh_interp(t(mask3),x1,x2,y1,y2,s1,s2,0,0);
end