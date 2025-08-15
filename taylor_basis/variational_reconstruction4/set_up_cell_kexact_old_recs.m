function CELLS = set_up_cell_kexact_old_recs(GRID,blk,idx_low,idx_high,degree,funs,use_bcs)
nvars = numel(funs);
IDX = generate_index_array(idx_low,idx_high);
CELLS = arrayfun( @(idx)set_up_cell(GRID,blk,idx,degree,nvars,funs,use_bcs), IDX );
% if (use_bcs)
%     error('need to implement this')
% end
end

function CELL = set_up_cell(GRID,blk,idx,degree,n_vars,funs1,use_bcs)

% q_scale = [1.0,75.0,75.0,75.0,100000.0];
% q_scale = [1,1,1,1,1];
% select = @(M,r) M(r);
% qfun1 = @(x)svf_prim(x,inputs);
% qfuns = {@(x)select(qfun1(x),1),@(x)select(qfun1(x),2),@(x)select(qfun1(x),3),@(x)select(qfun1(x),4),@(x)select(qfun1(x),5)};

funs = cell(numel(funs1),1);
for i = 1:numel(funs1)
    funs{i} = @(x) funs1{i}( x(1), x(2) );
end


[n_dim,~] = get_dim_and_nbor_info(GRID,blk,[idx{:}]);
% stencil_info
n_stencil = get_n_terms( n_dim, degree );
balanced   = true;
print_iter = false;

N = [GRID.gblock(1).imax,GRID.gblock(1).jmax,GRID.gblock(1).kmax];
block_info_list(1)           = block_info_t_old();
block_info_list(1).block_id  = 1;
block_info_list(1).Ncells(:) = max(N - 1,1);

% make stencil
% [ stencil, n_stencil, stencil_cells ] = build_stencil( block_id, idx, n_stencil, grid, block_info_list, balanced, print_iter )
[ stencil, ~, stencil_cells ] = cell_t_old.build_stencil( blk, [idx{:}], n_stencil, GRID, block_info_list, balanced, print_iter );

% soln_in_stencil = calc_soln_in_stencil_quads(n_var,funs,stencil,grid)
soln_stencil = calc_soln_in_stencil_quads(n_vars,n_dim,funs,stencil,GRID);

% set_bc_constraints( this, bound_num, n_var, n_dim, n_constraints, constraint_types, VAR_IDX, constraint_eval )
if (use_bcs)
    for b = 1:6
        % this = set_bc_constraints( this, bound_num, n_var, n_dim, n_constraints, constraint_types, VAR_IDX, constraint_eval )
        GRID.gblock(blk) = GRID.gblock(1).set_bc_constraints(b,n_vars,n_dim,n_vars,repmat(1000,[1,n_vars]),num2cell(1:n_vars),funs);
        % GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(1,5,2,1,3000,{[2,3,4]},[0]);
    end
    neighbor_max_degree = 0;
    mean_factor = 1000;
    boundary_factor = 1000;
    CELL = reconstruct_mod_t_old();
    CELL = CELL.setup_boundary( n_vars, n_dim, degree, stencil, GRID, neighbor_max_degree, stencil_cells, mean_factor, boundary_factor );
    CELL.polynomial = CELL.least_squares.solve(soln_stencil,CELL.polynomial,CELL.moments,stencil,GRID);
else
    CELL = reconstruct_mod_t_old();
    CELL = CELL.setup_interior( n_vars, n_dim, degree, stencil, GRID );
    CELL.polynomial = CELL.least_squares.solve( soln_stencil, CELL.polynomial, CELL.moments, stencil, GRID );
end
% pfun1 = @(x)rec1.polynomial.evaluate(x);


% pfun2 = @(x)rec2.polynomial.evaluate(x);

% CELL = reconstruct_mod_t_old();
% CELL = CELL.setup_boundary2( n_vars, n_dim, degree, stencil, GRID, neighbor_max_degree, stencil_cells );
% CELL.polynomial = CELL.least_squares.solve(soln_stencil,CELL.polynomial,CELL.moments,stencil,GRID);
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

function [face_quads,normals,face_nbor_idx,bc_idx] = get_face_info(GRID,blk,idx,n_dim)
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
bc_idx        = zeros(n_faces,1);

cnt = 0;
cnt_bc = 0;
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
    else
        cnt_bc = cnt_bc + 1;
        bc_idx(cnt_bc) = i;
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

bc_idx = bc_idx(1:cnt_bc);

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

function n_terms = get_n_terms( dimen, total_degree )
if (total_degree < 0)
    n_terms = 0;
else
    n_terms = nchoosek( dimen + total_degree, total_degree );
end
end

function soln_in_stencil = calc_soln_in_stencil_quads(n_var,n_dim,funs,stencil,grid)
n_stencil = size(stencil,2);
soln_in_stencil = zeros(n_stencil,n_var);
for si = 1:n_stencil
    b   = stencil(1,si);
    idx = stencil(2:4,si);
    quad = grid.gblock(b).grid_vars.quad(idx(1),idx(2),idx(3));
    vol  = grid.gblock(b).grid_vars.volume(idx(1),idx(2),idx(3));
    tmp_soln = zeros(n_var,quad.n_quad);
    for n = 1:quad.n_quad
        for vi = 1:n_var
            % tmp_pt = num2cell(quad.quad_pts(:,n));
            % tmp_soln(vi,n) = funs{vi}( tmp_pt{1:n_dim} );
            tmp_pt = quad.quad_pts(:,n);
            tmp_soln(vi,n) = funs{vi}( tmp_pt );
        end
    end
    soln_in_stencil(si,:) = quad.integrate(tmp_soln)/vol;
end
end