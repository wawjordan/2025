classdef rec_t
    properties
        n_blocks      (1,1) {mustBeInteger}
        n_dim         (1,1) {mustBeInteger}
        degree        (1,1) {mustBeInteger}
        n_vars        (1,1) {mustBeInteger}
        b             (:,1) rec_block_t
    end
    methods
        function this = rec_t( grid, n_dim, degree, n_vars, ext_fun )
            if ( nargin > 0 )
                this.n_dim  = n_dim;
                this.n_vars = n_vars;
                this.degree = degree;
                this.n_blocks = grid.nblocks;
                this.b = repmat(rec_block_t(),this.n_blocks);
                for i = 1:this.n_blocks
                    this.b(i) = rec_block_t( grid, i, n_dim, degree, n_vars );
                end

                if ( nargin == 5 )
                    for i = 1:this.n_blocks
                        this.b(i) = this.b(i).set_cell_avgs(grid.gblock(i),n_vars,(1:n_vars),ext_fun);
                    end
                end
            end
        end

        function err_norms = get_block_err_norms( this, grid, blk, var_idx, n_terms, norms, eval_fun )
            max_L_norm = 10;
            n_var  = numel(var_idx);
            n_norm = numel(norms);
            err_norms(1:n_var,1:n_norm) = 0.0;
            for n = 1:n_norm
                if (norm(n)>max_L_norm)
                    for i = 1:this.b(blk).n_cells_total
                        cell_idx = monomial_basis_t.global2local(i,this.b(blk).n_cells(1:this.n_dim));
                        tmp_idx(1:3) = 1;
                        tmp_idx(1:this.n_dim) = cell_idx;
                        quad = grid.gblock(blk).grid_vars.quad( tmp_idx(1), tmp_idx(2), tmp_idx(3) );
                        err_norms(:,n) = max(err_norms(:,n), this.b(blk).get_cell_error( quad, i, n_terms, norms(n), n_var, var_idx, eval_fun ));
                    end
                else
                    for i = 1:this.b(blk).n_cells_total
                        cell_idx = monomial_basis_t.global2local(i,this.b(blk).n_cells(1:this.n_dim));
                        tmp_idx(1:3) = 1;
                        tmp_idx(1:this.n_dim) = cell_idx;
                        quad = grid.gblock(blk).grid_vars.quad( tmp_idx(1), tmp_idx(2), tmp_idx(3) );
                        err_norms(:,n) = err_norms(:,n) + this.b(blk).get_cell_error( quad, i, n_terms, norms(n), n_var, var_idx, eval_fun );
                    end
                end
                err_norms(:,n) = err_norms(:,n) / this.b(blk).n_cells_total;
            end
        end

        function this = solve_rec(this, grid, ext_fun )
            var_idx = (1:this.n_vars);
            term_start = 1;
            for blk = 1:this.n_blocks
                term_end = this.b(blk).p.n_terms;
                this.b(blk) = this.b(blk).init_cells(grid,term_start,term_end);
                if ( nargin == 3)
                    this.b(blk) = this.b(blk).set_cell_avgs( grid.gblock(blk), this.n_vars, var_idx, ext_fun );
                end
                this.b(blk) = this.b(blk).solve_block( term_end, this.n_vars, var_idx );
                if ( nargin == 3)
                    error_norms = this.get_block_err_norms( grid, blk, var_idx, term_end, [1,2,realmax], ext_fun );
                    fprintf('Block %d error norms:\n',blk);
                    for v = 1:this.n_vars
                        fprintf('%d, %g, %g %g\n',v,error_norms(v,:));
                    end
                end
            end
        end
    end
end