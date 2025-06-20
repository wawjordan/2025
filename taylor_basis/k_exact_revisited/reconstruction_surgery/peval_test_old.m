%% peval test
clc; clear; close all;
rng('default')

n_dim = 3;
degree = 2;
n_var  = 1;

p = polynomial_mod_t_old(n_var, n_dim, degree);
p.coeff = rand(p.n_var,p.n_terms);
% p.coeff(:,5:end) = 0;
x = rand(p.n_dim,1);


pval = p.evaluate(x);
dpval1 = p.gradient_old(x)
dpval2 = p.gradient(x)

d2pval2 = p.hessian(x)
n_quad = ceil((degree+1)/2);
quad_ref = quad_t_old.create_quad_ref_3D(n_quad);
moments = compute_grid_moments_mod_old(p.exponents,p.x_ref,quad_ref)


