function [xb,yb,tb] = move_boundary_pt(xb0,xb1,yb0,yb1,xp,yp,tm1,tp1,yfun)

% estimate slope
del = 1e-6;

dy = (yb1 - yb0);
dy = ( (dy>=0) - (dy<0) )*max(abs(dy),del);

dx = (xb1 - xb0);
dx = ( (dx>=0) - (dx<0) )*max(abs(dx),del);

m = dy/dx;
xm = 1/m;
xb = (yp - yb0 + m*xb0 + xm*xp) / (m + xm);

yb_guess = -xm*(xb - xp) + yp;
fun = @(t) (yb_guess - yfun(t)).^2;
tb = fminbnd(fun,tm1,tp1);
yb = yfun(tb);
end