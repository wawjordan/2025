function val = determine_polygon_orientation(x,y)
val = 0;
if (numel(x)<3)
    return
end
val = sign(polygon_area(x,y));
end