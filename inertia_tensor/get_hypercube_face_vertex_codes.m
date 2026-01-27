function [v1,v2] = get_hypercube_face_vertex_codes(n_dim,dir,id)
% demonstration function to show selection of all vertices of a given face on a hypercube
% inputs:
%   n_dim - number of dimensions
%     dir - axis perpendicular to face
%      id - high or low face
%
% outputs:
%    v1 - vertices using loops
%    v2 - vertices using bitget

% there are 2^(n_dim-1) vertices for a n_dim-hypercube

% preallocate
v1 = zeros(2^(n_dim-1),n_dim);

% set the indices that don't change
v1(:,dir) = id;

% set the indices that change
for n = 1:2^(n_dim-1)
    tmp = n-1;
    i = 0;
    for d = 1:n_dim-1
        i = i+1;
        if (i == dir)
            i = i + 1;
        end
        v1(n,i) = mod(tmp,2);
        tmp = floor(tmp/2);
    end
end

% preallocate
v2 = zeros(2^(n_dim-1),n_dim);
% set the indices that don't change
v2(:,dir) = id;
% set the indices that change
for n = 1:2^(n_dim-1)
    bits = bitget(n-1,1:n_dim-1);
    v2(n,1:dir-1)     = bits(1:dir-1);
    v2(n,dir+1:n_dim) = bits(dir:n_dim-1);
end
end