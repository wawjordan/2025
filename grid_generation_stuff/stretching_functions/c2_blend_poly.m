function [xi,dxi,ddxi] = c2_blend_poly(t)
mask1 = t<=0;
mask2 = t>=1;
mask3 = (~mask1)&(~mask2);
xi   = zeros(size(t));
xi(mask2) = 1;
xi(mask3)   = t(mask3).^3.*(t(mask3).*(  6*t(mask3) -  15) + 10);

if (nargout>1)
    dxi  = zeros(size(t));
    dxi(mask3)  = t(mask3).^2.*(t(mask3).*( 30*t(mask3) -  60) + 30);
end

if (nargout>2)
    ddxi = zeros(size(t));
    ddxi(mask3) = t(mask3)   .*(t(mask3).*(120*t(mask3) - 180) + 60);
end

end