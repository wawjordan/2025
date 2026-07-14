function [map,der] = moretti_map(point,h,delta,index,factor)
% from Moretti (1976)
% function [map,der] = moretti_map(index,j,point,factor)
if (index==-1)
    delta = 1.0/delta;
end

p = point;
d1 = real(h);
d2 = imag(h);
ci = 1;

if (index==1)
    p = (point-h)/d1 + 1.0;
end

xm = real(p) - 1.0;
ym = imag(p);
xp = xm + 2.0;
y2 = ym.^2;
rm2 = xm.^2 + y2;
rp2 = xp.^2 + y2;
thp = atan(ym/xp);
if (xm==0) % goto 1
    thm = fsign(pi/2,ym);
    d   = ( (rm2/rp2)^(0.5*delta) )*exp(ci*delta*(thm-thp));
    map = (1.0 + d)/(1.0 - d);
    if (index==-1) % goto 4
        map = h + (map-1.0)*d1;
        der = d1*(point^2-1.0)/(delta*(map-h)*(map*complex(d1,-d2)));
        return
    end
    der = d1*delta*(map^2-1.0)/((point-h)*(point*complex(d1,-d2)));
    return
end
thm = atan(ym/xm);
if (xm>0) % goto 2
    d   = ( (rm2/rp2)^(0.5*delta) )*exp(ci*delta*(thm-thp));
    map = (1.0 + d)/(1.0 - d);
    if (index==-1) % goto 4
        map = h + (map-1.0)*d1;
        der = d1*(point^2-1.0)/(delta*(map-h)*(map*complex(d1,-d2)));
        return
    end
    der = d1*delta*(map^2-1.0)/((point-h)*(point*complex(d1,-d2)));
    return
end
thm = thm*fsign(pi,factor);
% goto 2
d   = ( (rm2/rp2)^(0.5*delta) )*exp(ci*delta*(thm-thp));
map = (1.0 + d)/(1.0 - d);
if (index==-1) % goto 4
    map = h + (map-1.0)*d1;
    der = d1*(point^2-1.0)/(delta*(map-h)*(map*complex(d1,-d2)));
    return
end
der = d1*delta*(map^2-1.0)/((point-h)*(point*complex(d1,-d2)));
return
end

function s = fsign(a,b)
if (a>=0)
    s = b;
else
    s = -b;
end
end

%                       delta
%  zeta - 1   / z - h \
%  -------- = |-------|                                              (9)
%  zeta + 1   \ z + h*/
%
% where:
% h     is the complex coordinate of the corner or edge in the z-plane
% h*    is its complex conjugate
% delta is the inverse of the external angle in multiples of pi
% In the zeta-plane, zeta=+/-1 are the images of z=h / z=h*
% but the contour is now smooth at zeta+/-1
% All other sharp angles on the contour remain unaltered in the zeta-plane
%
% The mapping defined by (9) is, to within a scale factor, the classic
% Karman-Trefftz airfoil-generating function.

% Note, first, that (9) can be decomposed into two parts:
%                delta
%      / z - h \
%  d = |-------|                                                    (10)
%      \ z + h*/
%
%         1 + d
%  zeta = -----                                                     (11)
%         1 - d
%
% In turn, (10) can be computed in this form,
%                delta
%      | z - h |      /                                    \
%  d = |-------|  exp | i delta[( arg(z - h) - arg(z + h*)]|        (12)
%      | z + h*|      \                                    /
%
% Now
%                    / y - beta  \
%  arg(z - h) =  atan|-----------| = theta_M                        (13)
%                    \ x - alpha /
%                    / y - beta  \
%  arg(z - h*) = atan|-----------| = theta_P
%                    \ x + alpha /
%
% if       h = alpha + i beta                                       (14)

% COMPLEX POINT,HINGE,CI,P,D,DER,H(10,10)
% COMMON/CMAP/PI,PI02,CI,H,POWER(10)
% DER IS DERIVATIVE  D(Z(J+1))/D(Z(J))

% function [map,der] = moretti_map(index,j,point,factor)
% hinge = h(j,j);
% pow   = power(J);
% if (index==-1)
%     pow = 1.0/pow;
% end
% 
% p = point;
% d1 = real(hinge);
% d2 = imag(hinge);
% 
% if (index==1)
%     p = (point-hinge)/d1 + 1.0;
% end
% 
% xm = real(p) - 1.0;
% ym = imag(p);
% xp = xm + 2.0;
% y2 = ym.^2;
% rm2 = xm.^2 + y2;
% rp2 = xp.^2 + y2;
% thp = atan(ym/xp);
% if (xm==0) % goto 1
%     thm = sign(pi/2,ym);
%     d   = ( (rm2/rp2)^(0.5*pow) )*exp(ci*pow*(thm-thp));
%     map = (1.0 + d)/(1.0 - d);
%     if (index==-1) % goto 4
%         map = hinge + (map-1.0)*d1;
%         der = d1*(point^2-1.0)/(pow*(map-hinge)*(map*complex(d1,-d2)));
%         return
%     end
%     der = d1*pow*(map^2-1.0)/((point-hinge)*(point*complex(d1,-d2)));
%     return
% end
% thm = atan(ym/xm);
% if (xm>0) % goto 2
%     d   = ( (rm2/rp2)^(0.5*pow) )*exp(ci*pow*(thm-thp));
%     map = (1.0 + d)/(1.0 - d);
%     if (index==-1) % goto 4
%         map = hinge + (map-1.0)*d1;
%         der = d1*(point^2-1.0)/(pow*(map-hinge)*(map*complex(d1,-d2)));
%         return
%     end
%     der = d1*pow*(map^2-1.0)/((point-hinge)*(point*complex(d1,-d2)));
%     return
% end
% thm = thm*sign(pi,factor);
% % goto 2
% d   = ( (rm2/rp2)^(0.5*pow) )*exp(ci*pow*(thm-thp));
% map = (1.0 + d)/(1.0 - d);
% if (index==-1) % goto 4
%     map = hinge + (map-1.0)*d1;
%     der = d1*(point^2-1.0)/(pow*(map-hinge)*(map*complex(d1,-d2)));
%     return
% end
% der = d1*pow*(map^2-1.0)/((point-hinge)*(point*complex(d1,-d2)));
% return
% end