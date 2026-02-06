%% Generate Polynomial Derivatives for Verifying Reconstruction Derivatives
clc; clear; close all;

syms X(3)

% test_fun_sym(x,y) = 1.0 + y + x + y.^2 + x.*y + x.^2;
% test_fun_sym(x,y) = -4.0 + 2.1*y + 0.02*x + 0.006*y.^2 -0.2*x.*y + x.^2;
% test_fun = matlabFunction(test_fun_sym);
% test_fun_grad = matlabFunction(gradient(test_fun_sym,[x,y]));
% test_fun_hess = matlabFunction(hessian(test_fun_sym,[x,y]));
