function [rho,u,v,p] = simple_inviscid_taylor_green_vortex_2D_vel(x,y,L,xo,yo)
mapx = @(x) (pi/L)*(x-xo) - pi/2;
mapy = @(y) (pi/L)*(y-yo) - pi/2;
u =  sin(mapx(x)) .* cos(mapy(y));
v = -cos(mapx(x)) .* sin(mapy(y));
% w = 0 *x.*y;
rho = 0 *x.*y + 1;
p = (rho/4).*( cos(2*mapx(x)) + cos(2*mapy(y)) );

% varargout{1} = rho;
% varargout{2} = u;
% varargout{3} = v;
% varargout{4} = p;
end

