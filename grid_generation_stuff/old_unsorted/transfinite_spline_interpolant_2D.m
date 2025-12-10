function F = transfinite_spline_interpolant_2D(XY1,XY2,XY3,XY4)

G = struct();
XY = {XY1,XY2,XY3,XY4};
for i = 1:4
    N = size(XY{i},1);
    t = linspace(-1,1,N);
    G(i).G = struct();
    G(i).G.sx = my_cubic_spline(t,XY{i}(:,1),1.0e99,1.0e99);
    G(i).G.sy = my_cubic_spline(t,XY{i}(:,2),1.0e99,1.0e99);
    G(i).G.x = @(t) G(i).G.sx.eval(t);
    G(i).G.y = @(t) G(i).G.sy.eval(t);
end

F = transfinite_interpolant(G(1).G,G(2).G,G(3).G,G(4).G);

end