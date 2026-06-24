function [x_sorted,y_sorted] = sort_points_counterclockwise(x,y,c)

if ( nargin == 3 )
    xc = c(1);
    yc = c(2);
else
    xc = mean(x);
    yc = mean(y);
end

angles = atan2(y-yc,x-xc);
[~,idx] = sort(angles);
x_sorted = x(idx);
y_sorted = y(idx);
end

% function points_sorted = sort_points_clockwise(points)
%     center = [mean(points(:,1)),mean(points(:,2))];
%     angles = atan2(points(:,2) - center(2), points(:,1) - center(1) );
%     [~,idx] = sort(angles);
%     points_sorted = points(idx,:);
% end