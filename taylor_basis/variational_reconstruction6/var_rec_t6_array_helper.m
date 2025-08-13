classdef var_rec_t6_array_helper

    properties
        CELLS   (:,:,:) var_rec_t6
        GRID      (1,1) grid_type
        GRID_VARS (1,1) grid_vars6
        degree  (1,1) {mustBeInteger} = 0
        n_vars  (1,1) {mustBeInteger} = 0
        dist_opt (1,:) {mustBeText}   = 'scalar_dist'
    end

    methods
        function this = var_rec_t6_array_helper( GRID, blk, idx_low, idx_high, degree, funs, n1, dist_opt )
            if nargin<1
                return
            end
            this.GRID = GRID;
            this.GRID_VARS = grid_vars6( GRID, blk );
            this.degree = degree;
            this.n_vars = numel(funs);
            this.dist_opt = dist_opt;
            IDX = generate_index_array(idx_low(1:this.GRID_VARS.n_dim),idx_high(1:this.GRID_VARS.n_dim));
            this.CELLS = arrayfun( @(idx)this.set_up_cell(idx{:}), IDX );
            this.CELLS = arrayfun(@(CELL)CELL.set_cell_avg( funs, this.GRID_VARS.get_cell_quad(GRID,CELL.idx_slf) ), this.CELLS );
            this = this.setup_all(n1);
            this = this.set_all_RHS();
        end

        function CELL = set_up_cell(this,idx)
            if (numel(idx)>1)
                idx_slf = this.GRID_VARS.local_to_global_bnd(idx,this.GRID_VARS.lo,this.GRID_VARS.hi);
            else
                idx_slf = idx;
            end

            n_int   = this.GRID_VARS.n_int(idx_slf);
            n_bnd   = this.GRID_VARS.n_bnd(idx_slf);
            idx_int   = this.GRID_VARS.nbor_idx(idx_slf,1:n_int);
            f_idx_int = 1:n_int;
            f_idx_bnd = n_int+1:n_bnd;
            h_ref     = this.GRID_VARS.get_cell_scale(this.GRID,idx);
            volume_quad = this.GRID_VARS.get_cell_quad(this.GRID,idx);
            basis = zero_mean_basis6(this.GRID_VARS.n_dim,this.degree,h_ref,volume_quad);

            CELL = var_rec_t6( idx_slf,      ...
                               idx_int,      ...
                               f_idx_int,    ...
                               f_idx_bnd,    ...
                               this.n_vars,       ...
                               basis,        ...
                               this.dist_opt );
        end

        function A = construct_full_LHS_matrix( this )
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:this.GRID_VARS.ncells
                nnz = nnz + (this.CELLS(i).n_int + 1)*n_terms_local^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:this.GRID_VARS.ncells
                row_start = (i-1)*n_terms_local + 1;
                col_start = row_start;

                nbors  = this.CELLS( this.CELLS(i).idx_int );

                f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
                           this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];

                Ai = this.CELLS(i).get_self_LHS(  nbors, [f_quad{:}] );
                Bij = this.CELLS(i).get_nbor_LHS( nbors, [f_quad{:}] );

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
            n_var  = this.CELLS(1).n_vars;
            n_dim   = this.CELLS(1).basis.n_dim;
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:this.GRID_VARS.ncells
                nnz = nnz + (n_dim*2+1)*( n_terms_local * n_var )^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:this.GRID_VARS.ncells
                row_start = (i-1)*n_terms_local*n_var + 1;
                col_start = row_start;
                nbors  = this.CELLS( this.CELLS(i).idx_int );

                f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
                           this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];

                Ai1  = this.CELLS(i).get_self_LHS( nbors, [f_quad{:}] );
                Bij1 = this.CELLS(i).get_nbor_LHS( nbors, [f_quad{:}] );
                Ai = kron(eye(n_var),Ai1);
                [row_idx,col_idx,vals,nz_cnt] = var_rec_t6_array_helper.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,Ai);
                for n = 1:this.CELLS(i).n_int
                    Bij = kron(eye(n_var),Bij1(:,:,n));
                    j = this.CELLS(i).idx_int(n);
                    col_start = (j-1)*n_terms_local*n_var + 1;
                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t6_array_helper.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-Bij);
                end
            end
            A = sparse(row_idx(1:nz_cnt),col_idx(1:nz_cnt),vals(1:nz_cnt));
        end

        function b = construct_full_RHS( this )
            n_var  = this.CELLS(1).n_vars;
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_cells = this.GRID_VARS.ncells;
            b = zeros( n_cells*n_terms_local, n_var );
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
                           this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];
                RHS = this.CELLS(i).get_RHS( nbors, [f_quad{:}] );
                b(row_start:row_end,:) = RHS;
            end
        end

        function b = construct_full_RHS_fully_coupled( this )
            n1      = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_var  = this.CELLS(1).n_vars;
            b = zeros( this.GRID_VARS.ncells*n_terms_local*n_var,1);
            for i = 1:this.GRID_VARS.ncells
                row_start = (i-1)*n_terms_local*n_var + 1;
                row_end   =   i  *n_terms_local*n_var;
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
                           this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];
                RHS = this.CELLS(i).get_RHS( nbors, [f_quad{:}] );
                b(row_start:row_end) = RHS(:);
            end
        end

        function this = setup_all(this,n1)
            for i = 1:this.GRID_VARS.ncells
                nbors  = this.CELLS( this.CELLS(i).idx_int );
                f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
                           this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];
                this.CELLS(i) = this.CELLS(i).setup_rec( n1, nbors, [f_quad{:}] );
            end
        end

        function this = set_all_RHS( this )
            for i = 1:this.GRID_VARS.ncells
                nbors = this.CELLS( this.CELLS(i).idx_int );
                f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
                           this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];
                this.CELLS(i).RHS = this.CELLS(i).get_RHS( nbors, [f_quad{:}] );
            end
        end

        % function this = perform_iterative_reconstruction_SOR(this,omega,n_iter)
        %     n_var = this.CELLS(1).n_vars;
        % 
        %     for iter = 1:n_iter
        %         for i = 1:this.GRID_VARS.ncells
        %             n1      = this.CELLS(i).ord0;
        %             n_terms = this.CELLS(i).basis.n_terms;
        %             update = zeros(n_terms-n1,n_var);
        %             for n = 1:this.CELLS(i).n_int
        %                 j = this.CELLS(i).idx_int(n);
        %                 uj = this.CELLS(j).coefs(n1+1:n_terms,:);
        %                 for v = 1:n_var
        %                     update(:,v) = update(:,v) + this.CELLS(i).AinvB(:,:,n) * uj(:,v);
        %                 end
        %             end
        %             nbors = this.CELLS( this.CELLS(i).idx_int );
        %             % f_quad = [ this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,true);...
        %             %        this.GRID_VARS.get_face_quads(this.GRID,this.CELLS(i).idx_slf,false)];
        %             % bi    = this.CELLS(i).get_RHS( nbors, [f_quad{:}] );
        %             % for v = 1:n_var
        %             %     update(:,v) = update(:,v) + this.CELLS(i).self_LHS_inv * bi(:,v);
        %             % end
        %             update = update + this.CELLS(i).get_RHS_update( nbors );
        %             % this.CELLS(i).coefs(n1+1:n_terms,:) = ...
        %             %     (1-omega) * this.CELLS(i).coefs(n1+1:n_terms,:) ...
        %             %     + omega  * update;
        % 
        %             tmp = this.CELLS(i).coefs(n1+1:n_terms,:);
        %             tmp = (1-omega) * tmp + omega  * update;
        %             this.CELLS(i).coefs(n1+1:n_terms,:) = tmp;
        %         end
        %     end
        % end

        function this = perform_iterative_reconstruction_SOR(this,omega,n_iter)
            n_var = this.CELLS(1).n_vars;
            coefs = cell( this.GRID_VARS.ncells );
            for i = 1:this.GRID_VARS.ncells
                coefs{i} = this.CELLS(i).coefs;
            end
            for iter = 1:n_iter
                for i = 1:this.GRID_VARS.ncells
                    n1      = this.CELLS(i).ord0;
                    n_terms = this.CELLS(i).basis.n_terms;
                    update = zeros(n_terms-n1,n_var);
                    for n = 1:this.CELLS(i).n_int
                        j = this.CELLS(i).idx_int(n);
                        uj = coefs{j}(n1+1:n_terms,:);
                        for v = 1:n_var
                            update(:,v) = update(:,v) + this.CELLS(i).AinvB(:,:,n) * uj(:,v);
                        end
                    end
                    nbors = this.CELLS( this.CELLS(i).idx_int );
                    update = update + this.CELLS(i).get_RHS_update( nbors );
                    tmp = coefs{i}(n1+1:n_terms,:);
                    coefs{i}(n1+1:n_terms,:) = (1-omega) * tmp + omega  * update;
                end
            end

            for i = 1:this.GRID_VARS.ncells
                this.CELLS(i).coefs = coefs{i};
            end
        end

        function [this,LHS,RHS,coefs] = perform_reconstruction( this )
            n1            = this.CELLS(1).ord0;
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_var        = this.CELLS(1).n_vars;
            n_cells       = this.GRID_VARS.ncells;

            LHS = this.construct_full_LHS_matrix();
            RHS = this.construct_full_RHS();
            coefs = zeros( n_terms_local*n_cells,n_var);
            for v = 1:n_var
                coefs(:,v) = LHS\RHS(:,v);
            end

            for i = 1:this.GRID_VARS.ncells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                n_terms   = this.CELLS(i).basis.n_terms;
                this.CELLS(i).coefs(n1+1:n_terms,:) = coefs(row_start:row_end,:);
            end
        end

        function [this,LHS,RHS,coefs] = perform_reconstruction_fully_coupled( this )
            n1            = this.CELLS(1).ord0; 
            n_terms_local = this.CELLS(1).basis.n_terms - n1;
            n_var         = this.CELLS(1).n_vars;
            n_cells       = this.GRID_VARS.ncells;

            % update RHS coefficients if needed
            this = this.set_all_RHS();

            RHS = this.construct_full_RHS_fully_coupled();
            LHS = this.construct_full_LHS_matrix_fully_coupled();

            coefs = zeros(n_terms_local*n_cells,n_var);
            coefs_calc = LHS\RHS;
            for i = 1:n_cells
                n_terms   = this.CELLS(i).basis.n_terms;
                for j = 1:n_var
                    row_start = (i-1)*n_var*n_terms_local ...
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