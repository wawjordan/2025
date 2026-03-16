classdef rec_block_t
    properties
        block_num     (1,1) {mustBeInteger}
        n_dim         (1,1) {mustBeInteger}
        degree        (1,1) {mustBeInteger}
        n_vars        (1,1) {mustBeInteger}
        n_cells       (:,1) {mustBeInteger}
        n_cells_total (1,1) {mustBeInteger}
        p             (1,1) monomial_basis_t
        cells         (:,1) rec_cell_t
    end

    methods
        function this = rec_block_t( grid, block_num, n_dim, degree, n_vars )
            this.block_num        = block_num;
            this.n_dim            = n_dim;
            this.degree           = degree;
            this.n_vars           = n_vars;
            this.n_cells(1:n_dim) = grid.gblock(block_num).n_cells(1:n_dim);
            this.n_cells_total    = prod( this.n_cells );
            this.p                = monomial_basis_t(n_dim,degree);
            this.cells = repmat(rec_cell_t,this.n_cells_total,1);

            idx_tmp(1:3,1) = 1;
            for i = 1:this.n_cells_total
                idx_tmp(1:this.n_dim,1) = monomial_basis_t.global2local( i, this.n_cells );
                nodes = rec_block_t.get_nodes(grid,block_num,idx_tmp,n_dim);
                h_ref = rec_block_t.get_cell_h_ref( nodes, n_dim );
                quad = grid.gblock(block_num).grid_vars.quad( idx_tmp(1), idx_tmp(2), idx_tmp(3) );
                [n_nbors, nbor_block, nbor_idx, nbor_degree, n_bnd, bnd_idx ] = get_nbors( this, grid, i );
                this.cells(i) = rec_cell_t( this.p, ...
                                            block_num,  ...
                                            i,          ...
                                            n_nbors,    ...
                                            n_bnd,      ...
                                            nbor_block, ...
                                            nbor_idx,   ...
                                            nbor_degree,...
                                            bnd_idx,    ...
                                            n_var,      ...
                                            quad,       ...
                                            h_ref );
            end
        end
        function [n_nbors, nbor_block, nbor_idx, nbor_degree, n_bnd, bnd_idx ] = get_nbors( this, grid, lin_idx )
            sz_in = ceil( 1.5 * this.p.n_terms );
            [sz_out,nbor_block,nbor_idx,nbor_degree,on_boundary] = grow_stencil_basic( this.n_dim, grid, this.block_num, lin_idx, sz_in );
            n_nbors = sz_out;
            n_bnd = sum(on_boundary);
            bnd_idx = nbor_idx(on_boundary);
        end
    end

    methods (Static)
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
                offset = monomial_basis_t.global2local(n,nSub) - 1;
                for i = 1:n_dim
                    tmp_idx = num2cell( idx(:) + offset(:) );
                    nodes(i,n) = GRID.gblock(blk).(node_var{i})( tmp_idx{:} );
                end
            end
        end
        function [sz_out,nbor_block,nbor_idx,nbor_degree,on_boundary] = grow_stencil_basic( n_dim, grid, block_num, lin_idx, sz_in )
            blk        = block_num;
            idx(1:3,1) = 1;
            Ncells = grid.gblock(blk).Ncells;
            idx(1:n_dim,1) = monomial_basis_t.global2local( lin_idx, Ncells );
            n_stencil  = sz_in;
            balanced   = true;
            print_iter = false;

            % block information (needed to inform stencil growth near boundaries
            block_info_list(1)           = block_info_t();
            block_info_list(1).block_id  = block_num;
            block_info_list(1).Ncells(:) = grid.gblock(blk).Ncells;
            [ stencil_idx, n_stencil_out, stencil_cells ] = cell_t.build_stencil( blk, idx, n_stencil, grid, block_info_list, balanced, print_iter );
            sz_out = n_stencil_out - 1;
            nbor_block  = zeros(sz_out,1);
            nbor_idx    = zeros(sz_out,1);
            nbor_degree = zeros(sz_out,1);
            on_boundary = false(sz_out,1);

            nbor_block(1:sz_out,1) = stencil_idx(1,2:end);
            for n = 2:sz_out+1
                nbor_idx(n) = monomial_basis_t.local2global( stencil_idx(2:n_dim+1,n), Ncells(1:n_dim) );
                nbor_degree(n) = stencil_cells(n).degree;
                on_boundary(n) = any( stencil_cells(n).bnd(1:2*n_dim) );
            end
        end
    end
end