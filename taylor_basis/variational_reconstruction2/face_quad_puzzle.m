%% figuring out how to get the right face quad (06/03/2025)
clc; clear; close all;

idx = [1,1,1];

n_dim = 2;

n_faces = 2*n_dim;

face_idxs = ones(numel(idx),n_faces);
for i = 1:n_faces
    face_idx =  fix((i-1)/2)+1;
    offset = mod(i+2,2);
    factor =  2*offset-1;
    face_idxs(1:n_dim,i)  = idx(1:n_dim);
    face_idxs(face_idx,i) = face_idxs(face_idx,i)+offset;
end

face_idxs