classdef var_rec_t6_array_helper

    properties
        CELLS (:,:,:) var_rec_t6
        ncells        {mustBeInteger} = 0
    end

    methods
        function this = var_rec_t6_array_helper( GRID )
            if nargin<1
                return
            end
        end

        function A = construct_full_LHS_matrix( this )
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:this.ncells
                nnz = nnz + (this.CELLS(i).n_nbors + 1)*n_terms_local^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:this.ncells
                row_start = (i-1)*n_terms_local + 1;
                col_start = row_start;

                nbors  = this.CELLS( this.CELLS(i).idx_int );

                Ai = this.CELLS(i).get_self_LHS(  nbors, f_quad );
                Bij = this.CELLS(i).get_nbor_LHS( nbors, f_quad );

                [row_idx,col_idx,vals,nz_cnt] = var_rec_t6_array_helper.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,Ai);
                
                for n = 1:this.CELLS(i).n_int

                    j = this.CELLS(i).idx_int(n);

                    col_start = (j-1)*n_terms_local + 1;

                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t6_array_helper.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-Bij(:,:,n));
                end
            end
            A = sparse(row_idx(1:nz_cnt),col_idx(1:nz_cnt),vals(1:nz_cnt));
        end


        function A = construct_full_LHS_matrix_fully_coupled( this )
            n_vars  = this.CELLS(1).n_vars;
            n_dim   = this.CELLS(1).basis.n_dim;
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:this.ncells
                nnz = nnz + (n_dim*2+1)*( n_terms_local * n_vars )^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:this.ncells
                row_start = (i-1)*n_terms_local*n_vars + 1;
                col_start = row_start;
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                Ai1  = this.CELLS(i).get_self_LHS( nbors, f_quad );
                Bij1 = this.CELLS(i).get_nbor_LHS( nbors, f_quad );
                Ai = kron(eye(n_vars),Ai1);
                [row_idx,col_idx,vals,nz_cnt] = var_rec_t6_array_helper.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,Ai);
                for n = 1:this.CELLS(i).n_int
                    Bij = kron(eye(n_vars),Bij1(:,:,n));
                    j = this.CELLS(i).idx_int(n);
                    col_start = (j-1)*n_terms_local*n_vars + 1;
                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t6_array_helper.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-Bij);
                end
            end
            A = sparse(row_idx(1:nz_cnt),col_idx(1:nz_cnt),vals(1:nz_cnt));
        end

        function b = construct_full_RHS( this )
            n_vars  = this.CELLS(1).n_vars;
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_cells = this.ncells;
            b = zeros( n_cells*n_terms_local, n_vars );
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                RHS = this.CELLS(i).get_RHS( nbors, f_quad );
                b(row_start:row_end,:) = RHS;
            end
        end

        function b = construct_full_RHS_fully_coupled( this )
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_vars  = this.CELLS(1).n_vars;
            b = zeros( this.ncells*n_terms_local*n_vars,1);
            for i = 1:this.ncells
                row_start = (i-1)*n_terms_local*n_vars + 1;
                row_end   =   i  *n_terms_local*n_vars;
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                RHS = this.CELLS(i).get_RHS( nbors, f_quad );
                b(row_start:row_end) = RHS(:);
            end
        end

        function this = setup_all(this,n1)
            for i = 1:this.ncells
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                this.CELLS(i) = this.CELLS(i).setup_rec( n1, nbors );
            end
        end

        function this = set_all_RHS( this )
            for i = 1:this.ncells
                nbors = this.CELLS( this.CELLS(i).idx_int );
                this.CELLS(i).RHS = this.CELLS(i).get_RHS( nbors, f_quad );
            end
        end

        function this = perform_iterative_reconstruction_SOR(this,omega,n_iter)
            n_vars = this.CELLS(1).n_vars;
            
            for iter = 1:n_iter
                for i = 1:this.ncells
                    n1      = this.CELLS(i).ord0;
                    n_terms = this.CELLS(i).basis.n_terms;
                    update = zeros(n_terms-n1,n_vars);
                    for n = 1:this.CELLS(i).n_int
                        j = this.CELLS(i).idx_int(n);
                        uj = this.CELLS(j).coefs(n1+1:n_terms,:);
                        for v = 1:n_vars
                            update(:,v) = update(:,v) + this.CELLS(i).AinvB(:,:,n) * uj(:,v);
                        end
                    end
                    nbors = this.CELLS( this.CELLS(i).idx_int );
                    bi    = this.CELLS(i).get_RHS( nbors, f_quad );
                    for v = 1:n_vars
                        update(:,v) = update(:,v) + this.CELLS(i).self_LHS_inv * bi(:,v);
                    end
                    this.CELLS(i).coefs(n1+1:n_terms,:) = ...
                        (1-omega) * this.CELLS(i).coefs(n1+1:n_terms,:) ...
                        + omega  * update;
                end
            end
        end

        function [this,LHS,RHS,coefs] = perform_reconstruction( this )
            n1            = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_vars        = this.CELLS(1).n_vars;
            n_cells       = this.ncells;

            LHS = this.construct_full_LHS_matrix();
            RHS = this.construct_full_RHS();
            coefs = zeros( n_terms_local*n_cells,n_vars);
            for v = 1:n_vars
                coefs(:,v) = LHS\RHS(:,v);
            end

            for i = 1:this.ncells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                n_terms   = this.CELLS(i).basis.n_terms;
                this.CELLS(i).coefs(n1+1:n_terms,:) = coefs(row_start:row_end,:);
            end
        end

        function [this,LHS,RHS,coefs] = perform_reconstruction_fully_coupled( this )
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_vars        = this.CELLS(1).n_vars;
            n_cells       = this.ncells;

            % update RHS coefficients if needed
            this.CELLS = this.set_all_RHS();

            RHS = this.construct_full_RHS_fully_coupled();
            LHS = var_rec_t6_array_helper.construct_full_LHS_matrix_fully_coupled();

            coefs = zeros(n_terms_local*n_cells,n_vars);
            coefs_calc = LHS\RHS;
            for i = 1:n_cells
                n_terms   = this.CELLS(i).basis.n_terms;
                for j = 1:n_vars
                    row_start = (i-1)*n_vars*n_terms_local ...
                              + (j-1)*n_terms_local + 1;
                    row_end   =  row_start + n_terms_local - 1;
                    this.CELLS(i).coefs(n1+1:n_terms,j) = coefs_calc(row_start:row_end);
                end
                % for debug & comparison
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                coefs(row_start:row_end,:) = this.CELLS(i).coefs(n1+1:n_terms,:);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function this = set_up_cell_var_recs( this, GRID, blk, idx_low, idx_high, degree,funs,n1,use_bcs,dist_opt)
        %     nvars = numel(funs);
        %     IDX = generate_index_array(idx_low,idx_high);
        % 
        %     this.CELLS = arrayfun( @(idx)set_up_cell(   GRID, ...
        %         blk, ...
        %         idx, ...
        %         degree, ...
        %         nvars, ...
        %         dist_opt),   IDX );
        %     CELLS = arrayfun(@(CELL)CELL.set_cell_avg(funs),CELLS);
        %     if (use_bcs)
        %         CELLS = var_rec_t6_array_helper.setup_all(n1,CELLS,funs);
        %     else
        %         CELLS = var_rec_t6_array_helper.setup_all(n1,CELLS,[]);
        %     end
        % end

        

        
    end

    methods (Static)
        function [row_idx,col_idx,vals,nz_cnt] = add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,array)
            [row_sz,col_sz] = size(array);
            cnt = nz_cnt;
            for j = 1:col_sz
                for i = 1:row_sz
                    cnt = cnt + 1;
                    row_idx(cnt) = row_start-1+i;
                    col_idx(cnt) = col_start-1+j;
                    vals(cnt)    = vals(cnt) + array(i,j);
                end
            end
            nz_cnt = cnt;
        end
        function val = in_bounds(idx,bnd_min,bnd_max)
            val = ( all(idx(:)>=bnd_min(:)) && all(idx(:)<=bnd_max(:)) || ...
                all(idx(:)<=bnd_min(:)) && all(idx(:)>=bnd_max(:)) );
        end
    end
end