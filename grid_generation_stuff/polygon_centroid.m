function c = polygon_centroid(x,y)
x1 = [x(:);x(1)];
y1 = [y(:);y(1)];
% A = 0.5*sum( x1(1:end-1).*y1(2:end) - x1(2:end).*y1(1:end-1) );
c(1,1) = sum( (x1(1:end-1) + x1(2:end) ) ...
            .*(x1(1:end-1).*y1(2:end) - x1(2:end).*y1(1:end-1) ) );
c(2,1) = sum( (y1(1:end-1) + y1(2:end) ) ...
            .*(x1(1:end-1).*y1(2:end) - x1(2:end).*y1(1:end-1) ) );
c = c / ( 3*sum( x1(1:end-1).*y1(2:end) - x1(2:end).*y1(1:end-1) ) );
end