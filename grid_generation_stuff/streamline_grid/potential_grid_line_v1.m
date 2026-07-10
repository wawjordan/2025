function [xfun,yfun] = potential_grid_line_v1(airfoil,x,y)
% find potential at location (x,y)
phi  = airfoil.phi_from_xy(x,y);
psi0 = airfoil.psi_from_xy(x,y);
xfun1  = @(y) airfoil.phi_x_from_y(phi,y,x);
psifun = @(y) airfoil.psi_from_xy(xfun1(y),y);
options = optimset('TolFun',1e-15,'TolX',1e-17);
obj_fun = @(psi,y) psifun(y) - psi;
yfun = @(psi) arrayfun(@(psi)fzero(@(y)obj_fun(psi,y),psi0,options),psi);
xfun = @(psi) xfun1(yfun(psi));
end

