function [x,y] = extrude_surface_pts(xs,ys,i1,i2,h)
s = determine_polygon_orientation(xs(i1:i2),ys(i1:i2));
if s == 0
    error('could not determine normal vector')
end
dx = gradient(xs(:));
dy = gradient(ys(:));
mag = sqrt( dx.^2 + dy.^2 );

% outward facing normal
n1 =  s*dy./mag;
n2 = -s*dx./mag;

if nargin<3
    h = mag;
    h1 = min(mag);
    h = h1*log((h/h1)*(exp(1)-1) + 1);
    % h = h1*sqrt(h/h1);
end
% add offset
x =xs + n1.*h(:);
y =ys + n2.*h(:);
end