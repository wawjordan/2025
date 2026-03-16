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
            if ( nargin > 0 )
                this.block_num        = block_num;
                this.n_dim            = n_dim;
                this.degree           = degree;
                this.n_vars           = n_vars;
                this.n_cells(1:n_dim) = grid.gblock(block_num).Ncells(1:n_dim);
                this.n_cells_total    = prod( this.n_cells );
                this.p                = monomial_basis_t(n_dim,degree);
                this.cells = repmat(rec_cell_t,this.n_cells_total,1);
    
                idx_tmp(1:3,1) = 1;
                for i = 1:this.n_cells_total
                    idx_tmp(1:this.n_dim,1) = monomial_basis_t.global2local( i, this.n_cells );
                    nodes = rec_block_t.get_nodes(grid,block_num,idx_tmp,n_dim);
                    h_ref = rec_block_t.get_max_cell_extents( nodes, n_dim );
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
                                                n_vars,     ...
                                                quad,       ...
                                                h_ref );
                end
            end
        end

        function [n_nbors, nbor_block, nbor_idx, nbor_degree, n_bnd, bnd_idx ] = get_nbors( this, grid, lin_idx )
            % sz_in = ceil( 1.5 * this.p.n_terms );
            sz_in = fix( 1.5 * this.p.n_terms );
            [sz_out,nbor_block,nbor_idx,nbor_degree,on_boundary] = rec_block_t.grow_stencil_basic( this.n_dim, grid, this.block_num, lin_idx, sz_in );
            n_nbors = sz_out;
            n_bnd = sum(on_boundary);
            bnd_idx = nbor_idx(on_boundary);
        end

        function [LHS,LHS_m,LHS_n] = get_cell_LHS(this,lin_idx,term_end,grid)
            i = lin_idx;
            LHS_m = this.cells(i).n_nbor;
            LHS_n = term_end-1;
            LHS = zeros(LHS_m,LHS_n);
            for j = 1:LHS_m
                shifted_moments = this.cells(i).basis.compute_shifted_moments( this.p, this.cells( this.cells(i).nbor_idx(j) ).basis );
                tmp_idx(1:3) = 1;
                tmp_idx(1:this.n_dim) = monomial_basis_t.global2local(this.cells(i).nbor_idx(j),this.n_cells);
                quad = grid.gblock( this.cells(i).nbor_block(j) ).grid_vars.quad( tmp_idx(1), tmp_idx(2), tmp_idx(3) );
                shifted_moments_q = this.cells(i).basis.compute_shifted_moments_quad( this.p, quad );
                LHS(j,:) = shifted_moments(2:term_end);
            end
        end

        function [RHS,n_nbors] = get_cell_RHS( this, lin_idx, n_vars, var_idx )
            i = lin_idx;
            n_nbors = this.cells(i).n_nbor;
            RHS = zeros(n_nbors,n_vars);
            for k = 1:n_nbors
                j = this.cells(i).nbor_idx(k);
                for v = 1:n_vars
                    RHS(k,var_idx(v)) = this.cells(j).coefs(1,var_idx(v)) - this.cells(i).coefs(1,var_idx(v));
                end
            end
        end

        function [LHS,Ainv] = init_cell(this,lin_idx,term_end,grid)
            i = lin_idx;
            [LHS,m,n] = this.get_cell_LHS(i,term_end,grid);
            Ainv(1:n,1:m) = pinv(LHS);
        end

        function this = init_cells(this,grid,~,term_end)
            for i = 1:this.n_cells_total
                [this.cells(i).LHS,this.cells(i).Ainv] = this.init_cell(i,term_end,grid);
            end
        end

        function coefs = solve_cell(this,lin_idx,term_end,n_var,var_idx)
            i = lin_idx;
            [RHS,n_nbors] = get_cell_RHS( this, i, n_var, var_idx );
            coefs = this.cells(i).coefs;
             for v = 1:n_var
                for n = 2:term_end
                    coefs(n,var_idx(v)) = dot( this.cells(i).Ainv(n-1,:),RHS(1:n_nbors,v) );
                end
             end
             tmp = this.cells(i).LHS * coefs(2:end);
             tmp2 = tmp - RHS;
        end

        function this = solve_block(this,term_end,n_var,var_idx)
            for i = 1:this.n_cells_total
                this.cells(i).coefs = this.solve_cell(i,term_end,n_var,var_idx);
            end
        end

        function this = set_cell_avgs(this,gblock,n_var,var_idx,eval_fun)
            tmp_idx(1:3) = 1;
            for i = 1:this.n_cells_total
                tmp_idx(1:this.n_dim) = monomial_basis_t.global2local(i,this.n_cells);
                quad = gblock.grid_vars.quad(tmp_idx(1),tmp_idx(2),tmp_idx(3));
                avg = rec_cell_t.get_cell_avg( quad, this.n_dim, n_var, eval_fun );
                for v = 1:n_var
                    this.cells(i).coefs(1,var_idx(v)) = avg(v);
                end
            end
        end

        function val = evaluate_block_rec( this, cell_idx, x, vars, n_terms, order )
            lin_idx = monomial_basis_t.local2global(cell_idx,this.n_cells);
            n_var  = numel(vars);

            if ( nargin == 6 )
                val = this.cells(lin_idx).evaluate_reconstruction_derivative( this.p, x, n_terms, n_var, vars, order );
            else
                val = this.cells(lin_idx).evaluate_reconstruction( this.p, x, n_terms, n_var, vars );
            end
        end

        function err = get_cell_error( this, quad, lin_idx, n_terms, norm, n_var, var_idx, eval_fun )
            max_L_norm = 10;
            tmp_val = zeros(n_var,quad.n_quad);

            for n = 1:quad.n_quad
                exact = eval_fun( quad.quad_pts(1:this.p.n_dim,n) );
                reconstructed = this.cells(lin_idx).evaluate_reconstruction( this.p, quad.quad_pts(:,n), n_terms, n_var, var_idx );
                tmp_val(:,n) = abs( reconstructed - exact );
            end
            if ( norm > max_L_norm )
                err = max(tmp_val,[],2);
            else
                err = ( quad.integrate( tmp_val.^norm ) ).^(1/norm);
                err = err / sum(quad.quad_wts)^(1/norm);
            end
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
            Ncells = grid.gblock(blk).Ncells(1:n_dim);
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
                nbor_idx(n-1) = monomial_basis_t.local2global( stencil_idx(2:n_dim+1,n), Ncells(1:n_dim) );
                nbor_degree(n-1) = stencil_cells(n).degree;
                on_boundary(n-1) = any( stencil_cells(n).bnd(1:2*n_dim) );
            end
        end

    end
end