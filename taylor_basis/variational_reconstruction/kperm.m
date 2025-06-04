function [p,p2] = kperm(n,k)
n_dim = numel(n);

diff = max(n(:) - k(:) + 1,0);
% if (any(diff<1))
%     p = 1;
%     return
% end
p = 1;
for d = 1:n_dim
    if k(d)>n(d)
        continue
    end
    for i = n(d):-1:n(d)-k(d)+1
        p = p * i;
    end
end

p2 = prod(factorial(n)) / prod( factorial(max(n-k,0)));

end

% fact_coef2 = prod(factorial(this.terms(n).exponents)) / prod( factorial( max(this.terms(n).exponents - order,0) ) );