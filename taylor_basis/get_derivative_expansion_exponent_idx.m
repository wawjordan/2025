function [terms,term_exp,exponents] = get_derivative_expansion_exponent_idx(n_dim,degree,term)
% terms:     [n_dim,1]
% exponents: [n_terms,n_dim+1]

terms = zeros(n_dim,1);
term_exp = zeros(n_dim+1,n_dim + 1);


[exp,idx,exponents] = get_exponents( n_dim, degree );

curr_exponent = exp(:,term);
term_degree   = sum(curr_exponent);

if (term_degree+1>degree)
    error('term must be less than terms with total degree')
end

start_idx = idx(term_degree+1);
stop_idx  = idx(term_degree+2);

term_exp(1,1)     = term;
term_exp(1,2:end) = curr_exponent;
cnt = 0;
for i = start_idx:stop_idx
    if (sum( abs(exp(:,i) - curr_exponent))==1)
        cnt = cnt + 1;
        terms(cnt) = i;
        term_exp(cnt+1,1) = terms(cnt);
        term_exp(cnt+1,2:end) = exp(:,i);
    end
end