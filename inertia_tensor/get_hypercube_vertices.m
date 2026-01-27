function [v1,v2] = get_hypercube_vertices(sz)
% demonstration function to enumerate all vertices of a hypercube
% inputs:
%   sz - size in each dimension
% outputs:
%    v1 - vertices using loops
%    v2 - vertices using bitget

n_dim = numel(sz);

% preallocate
v1 = zeros(2^n_dim,n_dim);
for n = 1:2^n_dim
    tmp = n-1;
    for i = 1:n_dim
        if ( mod(tmp,2)>0 )
            v1(n,i) = sz(i);
        end
        tmp = floor(tmp/2);
    end
end

% preallocate
v2 = zeros(2^n_dim,n_dim);
for n = 1:2^n_dim
    v2(n,:) = bitget(n-1,1:n_dim).*sz;
end

end