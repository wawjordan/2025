function [rho,u,v,p] = simple_inviscid_taylor_green_vortex_2D_vel_rot(x,y,L,xo,yo,theta)

% mapx = @(x,~)   x;
% mapy = @(~,y)   y;
mapx = @(x,y)   x * cos(theta) + y * sin(theta);
mapy = @(x,y)  -x * sin(theta) + y * cos(theta);


[rho,u1,v1,p] = simple_inviscid_taylor_green_vortex_2D_vel(mapx(x,y),mapy(x,y),L,xo,yo);

% u = u1;
% v = v1;
u =  u1 * cos(theta) - v1 * sin(theta);
v =  u1 * sin(theta) + v1 * cos(theta);

% varargout{1} = rho;
% varargout{2} = u;
% varargout{3} = v;
% varargout{4} = p;
end

