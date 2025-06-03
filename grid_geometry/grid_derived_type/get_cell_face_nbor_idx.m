function [nbor_idxs,cnt] = get_cell_face_nbor_idx(idx,bnd_min,bnd_max,max_dim)

% get number of dimensions
n_dim = numel(idx);

% number of "faces" for a given dimension
n_faces = 2*n_dim;

% allocate and initialize
nbor_idxs = ones(max_dim,n_faces+1);
for j = 1:n_faces
    nbor_idxs(1:n_dim,j) = idx;
end
cnt = 1;

% loop through dimensions, and add indices of valid face neighbors
for j = 1:n_dim
    for s = -1:2:1
        idx_tmp = idx;
        idx_tmp(j) = idx_tmp(j) + s;
        if ( in_bounds(idx_tmp,bnd_min,bnd_max) )
            cnt = cnt + 1;
            nbor_idxs(1:n_dim,cnt) = idx_tmp;
        end
    end
end
end

function val = in_bounds(idx,bnd_min,bnd_max)
val = ( all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:)) || ...
        all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:)) );
end