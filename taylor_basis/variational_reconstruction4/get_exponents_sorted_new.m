function E = get_exponents_sorted_new(n_dim,degree)

E = zeros(n_dim,nchoosek(n_dim+degree,degree));

% for n = 0:degree
%     for i = 1:nchoosek(n_dim+n,n)
%         for k = 1:n_dim
%             E(k,i) = 

if ( n_dim == 1 )
    for n = 0:degree
        E(n+1) = n;
    end
    return
end

tmp = zeros(n_dim,1);
full_degree = (degree+1)^n_dim;
p = full_degree;
for n = 0:full_degree
    tmp_deg = n + 1;
    for i = n_dim:-1:1
        p = fix( p/(degree+1) );
        deg = mod(tmp_deg-1,p)+1;
        tmp(i) = fix( (tmp_deg-deg)/p ) + 1;
        tmp_deg = deg;
        E(:,n+1) = tmp;
    end
end

end