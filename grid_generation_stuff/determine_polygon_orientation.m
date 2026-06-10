function val = determine_polygon_orientation(x,y)
N = numel(x);
val = 0;
if (N<3)
    return
end
% find vertex with smallest y (and largest x if tied)
[~,idx] = sortrows([y(:),x(:)],[1 2],["ascend" "descend"]);
% A = [x(idx(1));y(idx(1));0];
% B = [x(idx(2));y(idx(2));0];
% C = [x(idx(3));y(idx(3));0];
% val = sign(cross(B-A,C-A));
% sign(P) = sign( (x2-x1)(y3-y1) -(y2-y1)(x3-x1) )
val = sign( ( x(idx(2)) - x(idx(1)) ) ...
           *( y(idx(3)) - y(idx(1)) ) ...
           -( y(idx(2)) - y(idx(1)) ) ...
           *( x(idx(3)) - x(idx(1)) ) );
end