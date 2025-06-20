function monomial_eval = evaluate_monomial_d0_old( x, exponents, n_dim )
% eval (n_der == 0)
monomial_eval = 1.0;
monomial_eval = partial_monomial_evaluation_old(monomial_eval,x,exponents,n_dim);

end

% cnt = 0; for k = 1:3, for j = k:3, for i = j:3, cnt = cnt + 1; x = [0,0,0]; x(i) = x(i) + 1; x(j) = x(j) + 1; x(k) = x(k) + 1; fprintf('%2d | %d %d %d\n',cnt,x); end, end, end