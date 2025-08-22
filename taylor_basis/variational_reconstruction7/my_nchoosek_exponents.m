function [c,E] = my_nchoosek_exponents(n_dim,degree)

E = zeros( n_dim, nchoosek(n_dim+degree,degree) );

c = 0;
if (degree<0), return; end
c = 1;

cnt = 0;
for i = 1, min(n_dim,degree)
    c = c * ( n_dim + degree - (i-1) );
    c = c / i;
    for j = 1:c
        cnt = cnt + 1;
        E(cnt) = j;
    end
end
end