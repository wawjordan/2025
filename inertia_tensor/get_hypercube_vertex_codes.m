function [v1,v2] = get_hypercube_vertex_codes(n_dim)
% demonstration function to enumerate all vertices of a hypercube
% inputs:
%   n_dim - number of dimensions
% outputs:
%    v1 - vertices using loops
%    v2 - vertices using bitget

% preallocate
v1 = zeros(2^n_dim,n_dim);
for n = 1:2^n_dim
    tmp = n-1;
    for i = 1:n_dim
        v1(n,i) = mod(tmp,2);
        tmp = floor(tmp/2);
    end
end

% preallocate
v2 = zeros(2^n_dim,n_dim);
for n = 1:2^n_dim
    v2(n,:) = bitget(n-1,1:n_dim);
end

end

% function bits = binary_bits(d)
% bits = zeros(2^d,d);
% for n = 1:2^d
%     tmp = n-1;
%     for i = 1:d
%         bits(n,i) = mod(tmp,2);
%         tmp = floor(tmp/2);
%     end
% end
% end