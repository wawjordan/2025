function x = myquadspace( xmin, xmax, N )
% quadratic spacing ~> linear change in dx with increasing x
% but NOT linear change in density, which scales as 1/dx
x1 = linspace(0,1,N);
x = xmin + (x1.^2)*(xmax-xmin);
end