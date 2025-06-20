function monomial_eval = evaluate_monomial_d1_old( x, exponents, n_dim )
sz = n_dim;
monomial_eval = zeros(sz,1);
cnt = 0;
% gradient (n_der == 1)
for d = 1:n_dim
    cnt = cnt + 1;
    monomial_eval(cnt) = exponents(d);
    d_exponents = exponents;
    d_exponents(d) = d_exponents(d)-1;
    monomial_eval(cnt) = partial_monomial_evaluation_old(monomial_eval(cnt),x,d_exponents,n_dim);
end

end