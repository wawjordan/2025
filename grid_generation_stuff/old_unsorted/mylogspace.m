function x = mylogspace( xmin, xmax, beta, N )
% shifted logarithmic spacing ~> linear change in density with increasing x
% beta (> 1) is the base of the logarithm; increase for a greater change in grid
% density
x1 = linspace(beta,1,N);
x = xmax - log(x1)/log(beta)*(xmax-xmin); % change of base & scale
end