function monomial_eval = evaluate_monomial_d3( x, exponents, n_dim )
sz = nchoosek(n_dim+3-1,3);
monomial_eval = zeros(sz,1);
cnt = 0;
% 3rd derivative (only the symmetric part)
for d = 1:n_dim
    for di = d:n_dim
        for dj = di:n_dim
            cnt = cnt + 1;
            monomial_eval(cnt) = exponents(d);

            d_exponents = exponents;
    
            d_exponents(d)  = d_exponents(d)-1;
    
            monomial_eval(cnt) = monomial_eval(cnt)*d_exponents(di);
    
            d_exponents(di) = d_exponents(di)-1;

            monomial_eval(cnt) = monomial_eval(cnt)*d_exponents(di);
    
            d_exponents(dj) = d_exponents(dj)-1;
    
            monomial_eval(cnt) = partial_monomial_evaluation_old(monomial_eval(cnt),x,d_exponents,n_dim);
        end
    end
end

end
