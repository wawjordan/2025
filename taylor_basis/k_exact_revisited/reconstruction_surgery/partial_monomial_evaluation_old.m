function tmp_eval = partial_monomial_evaluation_old(tmp_eval,x,exponents,n_dim)
for d = 1:n_dim
    for i = exponents(d):-1:1
        tmp_eval = tmp_eval*x(d);
    end
end
end