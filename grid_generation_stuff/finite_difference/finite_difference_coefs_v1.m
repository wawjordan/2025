function [coefs,d_order] = finite_difference_coefs_v1(stencil,d_order)
stencil = stencil(:);

N = length(stencil);

if (d_order >= N)
    error('insufficient stencil points!')
    % coefs = [];
    % order = -1;
    % return
end



S = zeros(N,N);
for i = 1:N
    S(:,i) = stencil(i).^(0:N-1);
end

d = zeros(N,1);
d(d_order+1) = factorial(d_order);

coefs = S\d;

end