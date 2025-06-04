%% figuring out how to get the right face quad (06/03/2025)
clc; clear; close all;

idx = [1,1,1];
bnd_min = [1,1,1];
bnd_max = [5,5,1];

max_dim = numel(idx);
n_dim = 2;

[face_idxs,face_id] = get_cell_face_nbor_face_idx(idx,n_dim);
[nbor_idxs,cnt]     = get_cell_face_nbor_idx(idx,bnd_min,bnd_max,max_dim);

nbor_idxs = nbor_idxs(:,1:cnt);

% face_idxs - idx(:)

% for i = 1:n_faces
%     iG = local_to_global(face_idxs(i,:),)

function [face_idxs,face_id] = get_cell_face_nbor_face_idx(idx,n_dim)
n_faces = 2*n_dim;
face_idxs = ones(numel(idx),n_faces);
face_id   = 1:n_faces;
for i = 1:n_faces
    face_idx =  fix((i-1)/2)+1;
    offset = mod(i+2,2);
    % factor =  2*offset-1;
    face_idxs(1:n_dim,i)  = idx(1:n_dim);
    face_idxs(face_idx,i) = face_idxs(face_idx,i)+offset;
end
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