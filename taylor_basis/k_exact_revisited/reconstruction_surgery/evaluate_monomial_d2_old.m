function monomial_eval = evaluate_monomial_d2_old( x, exponents, n_dim )
sz = nchoosek(n_dim+2-1,2);
monomial_eval = zeros(sz,1);
cnt = 0;
% hessian (only the symmetric part)
for d = 1:n_dim
    for di = d:n_dim
        cnt = cnt + 1;
        monomial_eval(cnt) = exponents(d);

        d_exponents = exponents;

        d_exponents(d)  = d_exponents(d)-1;

        monomial_eval(cnt) = monomial_eval(cnt)*d_exponents(di);

        d_exponents(di) = d_exponents(di)-1;

        monomial_eval(cnt) = partial_monomial_evaluation_old(monomial_eval(cnt),x,d_exponents,n_dim);
    end
end

end