%% face indexing prototype 08/12/2025
clc; clear; close all;


bnd_min = [1,1,1];
bnd_max = [3,3,3];
max_dim = 3;

IDX = generate_index_array(bnd_min,bnd_max);

S(1).F = generate_index_array(bnd_min,bnd_max+[1,0,0]);
S(2).F = generate_index_array(bnd_min,bnd_max+[0,1,0]);
S(3).F = generate_index_array(bnd_min,bnd_max+[0,0,1]);


idx = [2,2,1];

% [nbor_idxs,cnt] = get_cell_face_nbor_idx(idx,bnd_min,bnd_max,max_dim);
% nbor_idxs(:,1:cnt)

linear_idx = local_to_global_bnd(idx,bnd_min,bnd_max);
[n_int,n_bnd,dir,face_idx] = get_cell_face_idx(linear_idx,bnd_min,bnd_max,max_dim);

for i = 1:n_int
    tmp = S( dir(i) ).F{abs(face_idx(i))};
    fprintf('%d : %d %d %d\n',i,tmp);
end
fprintf('\n')
for i = n_int+1:n_int+n_bnd
    tmp = S( dir(i) ).F{abs(face_idx(i))};
    fprintf('%d : %d %d %d\n',i,tmp);
end

function [n_int,n_bnd,dir,face_idx] = get_cell_face_idx(linear_idx,bnd_min,bnd_max,dim)
% get the subscripted cell index
cell_idx = global_to_local_bnd(linear_idx,bnd_min(1:dim),bnd_max(1:dim));

face_id = 1:2*dim;

[dir,linear_face_idx] = arrayfun(@(face_id)set_linear_face_idx( cell_idx, face_id, bnd_min(1:dim), bnd_max(1:dim) ),face_id);

interior = dir>0;

sort_idx = [face_id(interior),face_id(~interior)];

face_idx = linear_face_idx(sort_idx);
dir      = abs( dir(sort_idx) );
n_int = sum(interior);
n_bnd = 2*dim - n_int;
end

function [dir,offset,factor] = face_info(face_id)
dir    =  fix((face_id-1)/2)+1; % [1,1,2,2,3,3,...]
offset =  mod(face_id+1,2);     % [0,1]
factor =  2*offset-1;           % [-1, 1]
end

function [dir,linear_face_idx] = set_linear_face_idx(cell_idx,face_id,bnd_min,bnd_max)

[dir,offset,factor] = face_info(face_id);
idx_tmp = cell_idx;
idx_tmp(dir) = idx_tmp(dir) + offset;

nbor_idx = cell_idx;
nbor_idx(dir) = nbor_idx(dir)+factor;

interior = in_bounds(nbor_idx,bnd_min,bnd_max);

lo = bnd_min;
hi = bnd_max;
hi(dir) = hi(dir) + 1;
linear_face_idx = factor * local_to_global_bnd(idx_tmp,lo,hi);
dir(~interior) = -dir(~interior);
end

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

function iSub = global_to_local(iG,nSub)
nDims = numel(nSub);
iSub = zeros(1,nDims);
if (nDims==1)
    iSub(1) = iG;
    return
end
p = prod(nSub);
iGtmp = iG;
for i = nDims:-1:1
    p = fix( p/nSub(i) );
    iTmp = mod(iGtmp-1,p)+1;
    iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
    iGtmp = iTmp;
end
end

function iSub = global_to_local_bnd(iG,lo,hi)
% nSub = nSub + 2*nGhost;
% iSub = global_to_local(iG,nSub);
% iSub = iSub - nGhost;
nSub = hi(:) - lo(:) + 1;
iSub = global_to_local(iG,nSub);
iSub = iSub(:) + lo(:) - 1;

end

function iG = local_to_global(iSub,nSub)
iSub = iSub(:);
nSub = nSub(:);
nDims = numel(iSub);
p = 1;
iG = 1;
for i = 1:nDims
    iG = iG + ( iSub(i) - 1 )*p;
    p = p*nSub(i);
end

end

function iG = local_to_global_bnd(iSub,lo,hi)
% iSub = iSub(:) + nGhost(:);
iSub = iSub(:) - lo(:) + 1;
% nSub = nSub(:) + 2*nGhost(:);
nSub = hi(:) - lo(:) + 1;
iG = local_to_global(iSub,nSub);
end