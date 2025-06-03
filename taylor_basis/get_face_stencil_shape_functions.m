function [T,Q,id] = get_face_stencil_shape_functions(GRID,blk,idx,degree)
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
bnd_max = tmp_array(1:n_dim);
[nbor_idxs,cnt] = get_cell_face_nbor_idx(idx(1:n_dim),bnd_min,bnd_max,max_dim);

for i = 1:cnt
    nodes  = get_nodes(GRID,blk,nbor_idxs(:,i),n_dim);
    tmp_idx = num2cell( nbor_idxs(:,i) ); 
    Q(i,1) = GRID.gblock(blk).grid_vars.quad( tmp_idx{:} );
    T(i,1) = taylor_shape_functions(n_dim,degree,nodes,Q(i,1));
end

id = nbor_idxs(:,1:cnt);

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