function moments = compute_grid_moments_mod_old(exponents,x_ref,quad)

n_dim = size(exponents,1);
n_terms = size(exponents,2);
n_quad  = quad.n_quad;

moments = zeros(n_terms,1);

for m = 1:n_terms
    moment = zeros(1,n_quad);
    for n = 1:n_quad
        x = quad.quad_pts(:,n) - x_ref;
        moment(n) = evaluate_monomial_d0_old(x,exponents(:,m),n_dim);
    end
    moments(m) = quad.integrate(moment);
end

moments = moments / moments(1);

end