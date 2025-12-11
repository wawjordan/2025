function [xc,yc,A] = quad_grid_centroids(x,y)
% order = [1,2,4,2];
% x1 = x(1:end-1,1:end-1);
% x2 = x(2:end  ,1:end-1);
% x3 = x(2:end  ,2:end  );
% x4 = x(1:end-1,2:end  );
% y1 = y(1:end-1,1:end-1);
% y2 = y(2:end  ,1:end-1);
% y3 = y(2:end  ,2:end  );
% y4 = y(1:end-1,2:end  );

[xc,yc,A] = arrayfun( @(x1,x2,x3,x4,y1,y2,y3,y4)               ...
                        quad_centroid(x1,x2,x3,x4,y1,y2,y3,y4),...
                        x(1:end-1,1:end-1), ...
                        x(2:end  ,1:end-1), ...
                        x(2:end  ,2:end  ), ...
                        x(1:end-1,2:end  ), ...
                        y(1:end-1,1:end-1), ...
                        y(2:end  ,1:end-1), ...
                        y(2:end  ,2:end  ), ...
                        y(1:end-1,2:end  ) );
end

function [xc,yc,A] = quad_centroid(x1,x2,x3,x4,y1,y2,y3,y4)
[c,A] = poly_centroid([x1,x2,x3,x4],[y1,y2,y3,y4]);
xc = c(1);
yc = c(2);
end

function [c,A] = poly_centroid(x_coords,y_coords)
N = numel(x_coords);
c = [0;0];
A = 0;
for i = 1:N-1
    det = x_coords(i)*y_coords(i+1) - x_coords(i+1)*y_coords(i);
    A    = A    + det;
    c(1) = c(1) + ( x_coords(i) + x_coords(i+1) )*det;
    c(2) = c(2) + ( y_coords(i) + y_coords(i+1) )*det;
end
c = c / (6*abs(A));
end
% function A = poly_area(x_coords,y_coords)
% N = numel(x_coords);
% A = 0;
% for i = 1:N-1
%     A = A + x_coords(i)*y_coords(i+1) - x_coords(i+1)*y_coords(i);
% end
% A = 0.5 * A;
% end