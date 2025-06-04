function CELLS = set_up_cell_var_recs2(GRID,blk,idx_low,idx_high,degree,funs,n1)
nvars = numel(funs);
IDX = generate_index_array(idx_low,idx_high);
CELLS = arrayfun( @(idx)set_up_cell(GRID,blk,idx,degree,nvars), IDX );
CELLS = arrayfun(@(CELL)CELL.set_cell_avg(funs),CELLS);
CELLS = var_rec_t2.setup_all(n1,CELLS);
end

function CELL = set_up_cell(GRID,blk,idx,degree,n_vars)
[n_dim,nbor_idx] = get_dim_and_nbor_info(GRID,blk,[idx{:}]);

n_nbor = size(nbor_idx,2) - 1;
tmp_idx = num2cell( nbor_idx(:,1) );
Q = GRID.gblock(blk).grid_vars.quad( tmp_idx{:} );
[FQ,FN,nbor_face_id] = get_face_info(GRID,blk,idx,n_dim);
nodes = get_nodes(GRID,blk,nbor_idx(:,1),n_dim);
h_ref = get_max_cell_extents(nodes,n_dim);
T = zero_mean_basis2(n_dim,degree,h_ref,Q);
CELL = var_rec_t2(nbor_idx,nbor_face_id(1:n_nbor),T,Q,FQ,FN,n_vars);
end

function [n_dim,id] = get_dim_and_nbor_info(GRID,blk,idx)
max_dim = 3;
bnd_vars={'imax','jmax','kmax'};
tmp_array = zeros(max_dim,1);
n_dim = 0;
for i = 1:max_dim
    if isfield( GRID.gblock(blk), bnd_vars{i} ) ...
     || isprop( GRID.gblock(blk), bnd_vars{i} )
        if ( GRID.gblock(blk).(bnd_vars{i}) > 2 )
            n_dim = n_dim + 1;
            tmp_array(n_dim) = GRID.gblock(blk).(bnd_vars{i});
        end
    end
end

bnd_min = ones(n_dim,1);
bnd_max = tmp_array(1:n_dim)-1;
[nbor_idxs,cnt] = get_cell_face_nbor_idx(idx(1:n_dim),bnd_min,bnd_max,max_dim);
id = nbor_idxs(:,1:cnt);
end

function h_ref = get_max_cell_extents(node_list,n_dim)
h_min = min(node_list(1:n_dim,:),[],2);
h_max = max(node_list(1:n_dim,:),[],2);
h_ref = 0.5 * (h_max(:) - h_min(:));
end

function nodes = get_nodes(GRID,blk,idx,n_dim)
node_var={'x','y','z'};
max_dim = 3;
if (n_dim > max_dim)
    error('GRID type not set up for more than 3 dimensions')
end
n_nodes = 2^n_dim;
nodes = zeros(max_dim,n_nodes);
nSub = 2*ones(max_dim,1);
for n = 1:n_nodes
    offset = global_to_local(n,nSub) - 1;
    for i = 1:n_dim
        tmp_idx = num2cell( idx(:) + offset(:) );
        nodes(i,n) = GRID.gblock(blk).(node_var{i})( tmp_idx{:} );
    end
end
end

function [face_quads,normals,face_nbor_idx] = get_face_info(GRID,blk,idx,n_dim)
quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
norm_var={'xi_n','eta_n','zeta_n'};
max_dim = 3;
if (n_dim > max_dim)
    error('GRID type not set up for more than 3 dimensions')
end


n_faces = 2*n_dim;
face_quads = cell(n_faces,1);
normals    = cell(n_faces,1);

face_nbor_idx = zeros(n_faces,1);

cnt = 0;
for i = 1:n_faces
    face_idx =  fix((i-1)/2)+1; % [1,1,2,2,3,3]
    offset   =  mod(i+1,2);     % [0,1]
    factor   =  2*offset-1;     % [-1, 1]
    tmp_idx  =  [idx{:}];
    tmp_idx(face_idx) = tmp_idx(face_idx)+offset;

    nbor_idx = [idx{:}];
    nbor_idx(face_idx) = nbor_idx(face_idx)+factor;

    if ( in_bounds(nbor_idx(1:n_dim),ones(1,n_dim),GRID.gblock(blk).Ncells(1:n_dim)) )
        cnt = cnt + 1;
        face_nbor_idx(cnt) = i;
    end

    face_quads{i} = GRID.gblock(blk).grid_vars.(quad_var{face_idx}) ...
                                               ( tmp_idx(1),...
                                                 tmp_idx(2),...
                                                 tmp_idx(3) );
    normals{i}    = GRID.gblock(blk).grid_vars.(norm_var{face_idx}) ...
                                               ( tmp_idx(1),...
                                                 tmp_idx(2),...
                                                 tmp_idx(3) );
    normals{i}.v = normals{i}.v * factor;
end

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

function val = in_bounds(idx,bnd_min,bnd_max)
val = ( all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:)) || ...
        all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:)) );
end