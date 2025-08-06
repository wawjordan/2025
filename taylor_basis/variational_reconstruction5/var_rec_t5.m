classdef var_rec_t5

    properties
        n_nbors      (1,1)   {mustBeInteger} = 1 % number of face nbors
        idx          (:,1)   {mustBeInteger}     % self index
        nbor_idx     (:,:)   {mustBeInteger}     % nbor indices
        nbor_face_id (:,1)   {mustBeInteger}     % corresponding face ids
        bc_face_id   (:,1)   {mustBeInteger}     % face ids for bcs
        basis        (1,1)   zero_mean_basis5    % local basis functions
        quad         (1,1)   quad_t              % volume quadrature
        fquad        (:,1)   quad_t              % face quadrature
        fvec         (:,1)   face_vec_t          % face normal vecs
        n_vars       (1,1)   {mustBeInteger}     % number of eqns
        avgs         (:,1)   {mustBeNumeric}     % cell avgs
        coefs        (:,:)   {mustBeNumeric}     % polynomial coeffs
        self_LHS     (:,:)   {mustBeNumeric}     % LHS for self
        nbor_LHS     (:,:,:) {mustBeNumeric}     % LHS for each nbor
        RHS          (:,:)   {mustBeNumeric}     % RHS for self
        self_LHS_inv (:,:)   {mustBeNumeric}     % pseudo-inverse of self LHS
        AinvB        (:,:,:) {mustBeNumeric}     % product with self_LHS_inv
        Ainvb        (:,:)   {mustBeNumeric}
        c            (:,:)   {mustBeNumeric}     % RHS = c*u_bar + d
        d            (:,:)   {mustBeNumeric}
        Ainvc        (:,:)   {mustBeNumeric}
        Ainvd        (:,:)   {mustBeNumeric}
        get_nbor_dist
    end

    methods
        function this = var_rec_t5( nbor_idxs, ...
                                    nbor_face_id, ...
                                    bc_face_id, ...
                                    basis_funs, ...
                                    volume_quad, ...
                                    face_quads, ...
                                    face_normals, ...
                                    n_vars, ...
                                    dist_opt )
            if nargin<1
                return
            end
            this.idx      = nbor_idxs(:,1);
            this.nbor_idx = nbor_idxs(:,2:end);
            this.nbor_face_id = nbor_face_id;
            this.bc_face_id   = bc_face_id;
            this.n_nbors  = size(nbor_idxs,2) - 1;
            this.basis    = basis_funs;
            this.quad     = volume_quad;
            this.fquad    = [face_quads{:}];
            this.fvec     = [face_normals{:}];
            this.n_vars   = n_vars;
            this.avgs     = zeros(n_vars,1);
            this.coefs    = zeros(this.basis.n_terms,n_vars);
            this.self_LHS = zeros(this.basis.n_terms,this.basis.n_terms);
            this.self_LHS_inv = zeros(this.basis.n_terms,this.basis.n_terms);
            this.nbor_LHS = zeros(this.basis.n_terms,this.basis.n_terms,this.n_nbors);
            this.RHS      = zeros(this.basis.n_terms,this.n_vars);
            this.AinvB    = zeros(this.basis.n_terms,this.basis.n_terms,this.n_nbors);
            switch dist_opt
                case 'scalar_dist'
                    this.get_nbor_dist = @this.get_nbor_distance_scalar;
                case 'vector_dist'
                    this.get_nbor_dist = @this.get_nbor_distance_vec1;
                otherwise
                    error("Options are 'scalar_dist' or 'vector_dist'");
            end
        end

        function varargout = eval_reconstruction(this,x,n_terms)
            varargout = cell(this.n_vars,1);
            Bx = zeros(this.basis.n_terms,1);
            if nargin > 2
                for i = 1:n_terms
                    Bx(i) = this.basis.eval(i,x);
                end
            else
                for i = 1:this.basis.n_terms
                    Bx(i) = this.basis.eval(i,x);
                end
            end

            for v = 1:this.n_vars
                varargout{v} = sum(Bx.*this.coefs(:,v));
            end
        end

        function this = set_cell_avg(this,fun)
            n_dim  = this.basis.n_dim;
            n_quad = this.quad.n_quad;
            volume = this.quad.integrate(ones(1,n_quad));
            tmp_val = zeros(1,n_quad);
            for v = 1:numel(fun)
                for n = 1:n_quad
                   tmp_x      = num2cell(this.quad.quad_pts(1:n_dim,n));
                   tmp_val(n) = fun{v}(tmp_x{:});
                end
                this.avgs(v) = this.quad.integrate(tmp_val)/volume;
            end
            this.coefs(1,:) = this.avgs;
        end

        function this = setup_rec(this,n1,nbors,bc_funs)
            this.self_LHS = this.get_self_LHS(n1,nbors);
            if (~isempty(bc_funs))
            end
            this.self_LHS_inv = pinv(this.self_LHS);
            this.nbor_LHS = this.get_nbor_LHS(n1,nbors);
            this.RHS      = this.get_RHS(n1,nbors);
            n_nbor  = this.n_nbors;
            n_terms_local = this.basis.n_terms-n1;
            this.AinvB = zeros(n_terms_local,n_terms_local,n_nbor);
            for n = 1:n_nbor
                this.AinvB(:,:,n) = this.self_LHS_inv * this.nbor_LHS(:,:,n);
            end
        end

        function d = get_nbor_distance_scalar(this,nbor)
            dij = this.basis.x_ref - nbor.basis.x_ref;
            d   = sqrt( sum(dij.^2));
        end

        function d = get_nbor_distance_vec1(this,nbor)
            % d = (this.basis.h_ref + nbor.basis.h_ref)/2;
            d = abs(this.basis.x_ref - nbor.basis.x_ref);
        end

        function A = get_self_LHS(this,n1,nbors)
            n_nbor  = this.n_nbors;
            n_terms_local = this.basis.n_terms-n1;
            max_n_quad = 0;
            for k = 1:n_nbor
                max_n_quad = max(max_n_quad,this.fquad( this.nbor_face_id(k) ).n_quad );
            end
            d_basis = zeros(max_n_quad,n_terms_local+n1,n_terms_local);
            A = zeros(n_terms_local,n_terms_local);
            
            for k = 1:n_nbor
                dij_mag = this.get_nbor_distance_scalar(nbors(k));
                dij = this.get_nbor_dist(nbors(k));
                n_quad = this.fquad( this.nbor_face_id(k) ).n_quad;

                for q = 1:n_quad
                    point = this.fquad( this.nbor_face_id(k) ).quad_pts(:,q);
                    d_basis(q,:,:) = this.basis.calc_basis_derivatives(n1,point,dij);
                end
                A1 = zeros(n_terms_local,n_terms_local);
                for l = 1:n_terms_local
                    for m = 1:n_terms_local
                        sum_tmp = sum(d_basis(1:n_quad,:,l) .* d_basis(1:n_quad,:,m),2);
                        A1(m,l) = this.fquad( this.nbor_face_id(k) ).integrate( sum_tmp.' );
                    end
                end
                A = A + A1/dij_mag;
            end
        end

        function B = get_nbor_LHS(this,n1,nbors)
            n_nbor  = this.n_nbors;
            n_terms_local = this.basis.n_terms - n1;
            max_n_quad = 0;
            for k = 1:n_nbor
                max_n_quad = max(max_n_quad,this.fquad( this.nbor_face_id(k) ).n_quad );
            end
            d_basis_i = zeros(max_n_quad,n_terms_local+n1,n_terms_local);
            d_basis_j = zeros(max_n_quad,n_terms_local+n1,n_terms_local);

            B = zeros(n_terms_local,n_terms_local,n_nbor);

            for k = 1:n_nbor
                dij     = this.get_nbor_dist(nbors(k));
                dij_mag = this.get_nbor_distance_scalar(nbors(k));
                n_quad = this.fquad( this.nbor_face_id(k) ).n_quad;
                for q = 1:n_quad
                    point = this.fquad( this.nbor_face_id(k) ).quad_pts(:,q);
                    d_basis_i(q,:,:) = this.basis.calc_basis_derivatives(n1,point,dij);
                    d_basis_j(q,:,:) = nbors(k).basis.calc_basis_derivatives(n1,point,dij);
                end
                for r = 1:n_terms_local
                    for l = 1:n_terms_local
                        % sum of product of derivatives
                        sum_tmp = sum( d_basis_i(1:n_quad,:,l) .* d_basis_j(1:n_quad,:,r), 2 );
                        B(l,r,k) = this.fquad( this.nbor_face_id(k) ).integrate( sum_tmp.' );
                    end
                end
                B(:,:,k) = B(:,:,k) / dij_mag;
            end
        end

        function b = get_RHS(this,n1,nbors)
            n_nbor  = this.n_nbors;
            n_terms_local = this.basis.n_terms - n1;
            n_quad  = this.fquad( this.nbor_face_id(1) ).n_quad;

            b = zeros(n_terms_local,this.n_vars);

            for k = 1:n_nbor
                dij = this.get_nbor_dist(nbors(k));
                dij_mag = this.get_nbor_distance_scalar(nbors(k));
                L = zero_mean_basis5.get_length_scale(zeros(this.basis.n_dim,1),dij);
                avg_diff = nbors(k).avgs(:) - this.avgs(:);
                for l = 1:n_terms_local
                    int_tmp = zeros(1,n_quad);
                    for q = 1:n_quad
                        point = this.fquad( this.nbor_face_id(k) ).quad_pts(:,q);
                        int_tmp(q) = this.basis.eval( l+n1, point );
                    end
                    b(l,:) = b(l,:) + avg_diff(:).' * (1/dij_mag) * (L^2) * this.fquad( this.nbor_face_id(k) ).integrate(int_tmp);
                end
            end
        end


    end
    methods (Static)
        function A = construct_full_LHS_matrix(n1,CELLS)

            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            n_terms_local = CELLS(1).basis.n_terms - n1;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:n_cells
                nnz = nnz + (CELLS(i).n_nbors + 1)*n_terms_local^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local + 1;
                col_start = row_start;
                Ai = CELLS(i).self_LHS;
                [row_idx,col_idx,vals,nz_cnt] = var_rec_t5.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,Ai);
                n_nbors = CELLS(i).n_nbors;
                for n = 1:n_nbors
                    Bij = CELLS(i).nbor_LHS(:,:,n);

                    nbor_idx = CELLS(i).nbor_idx(1:n_dim,n);

                    j = zero_mean_basis5.local_to_global(nbor_idx,sz);

                    col_start = (j-1)*n_terms_local + 1;
                    
                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t5.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-Bij);
                end
            end
            A = sparse(row_idx(1:nz_cnt),col_idx(1:nz_cnt),vals(1:nz_cnt));
        end


        function A = construct_full_LHS_matrix_fully_coupled(n1,CELLS,var_idx)

            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            n_vars  = CELLS(1).n_vars;
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            n_terms_local = CELLS(1).basis.n_terms - n1;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:n_cells
                % nnz = nnz + (CELLS(i).n_nbors + 1)*( n_terms_local * n_vars )^2;
                nnz = nnz + (n_dim*2+1)*( n_terms_local * n_vars )^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local*n_vars + 1;
                col_start = row_start;
                Ai = kron(eye(n_vars),CELLS(i).self_LHS);
                [row_idx,col_idx,vals,nz_cnt] = var_rec_t5.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,Ai);
                n_nbors = CELLS(i).n_nbors;
                for n = 1:n_nbors
                    Bij = kron(eye(n_vars),CELLS(i).nbor_LHS(:,:,n));

                    nbor_idx = CELLS(i).nbor_idx(1:n_dim,n);

                    j = zero_mean_basis5.local_to_global(nbor_idx,sz);

                    col_start = (j-1)*n_terms_local*n_vars + 1;
                    
                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t5.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-Bij);
                end
                if (~isempty(var_idx))
                    % for n = 1:numel(CELLS(i).bc_face_id)
                    %     bc_fquad = CELLS(i).fquad( CELLS(i).bc_face_id(n) );
                    %     bc_fvec  = CELLS(i).fvec(  CELLS(i).bc_face_id(n) );
                    %     vrswt = var_rec_slip_wall_t(n1,CELLS(i),bc_fquad,bc_fvec,var_idx);
                    %     % for each sub-matrix
                    %     for j = 1:numel(var_idx)
                    %         row_start = (i-1)*n_vars*n_terms_local + (var_idx(j)-1)*n_terms_local + 1;
                    %         sub_row_start = (j-1)*n_terms_local + 1;
                    %         for k = 1:numel(var_idx)
                    %             col_start = (i-1)*n_vars*n_terms_local + (var_idx(k)-1)*n_terms_local + 1;
                    %             sub_col_start   = (k-1)*n_terms_local + 1;
                    %             sub_matrix = vrswt.E(sub_row_start:sub_row_start+n_terms_local-1,sub_col_start:sub_col_start+n_terms_local-1);
                    %             [row_idx,col_idx,vals,nz_cnt] = var_rec_t5.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,sub_matrix);
                    %         end
                    %     end
                    % end
                end
                
            end
            A = sparse(row_idx(1:nz_cnt),col_idx(1:nz_cnt),vals(1:nz_cnt));
        end

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

        function b = construct_full_RHS(n1,CELLS)
            n_cells = numel(CELLS);
            n_terms_local = CELLS(1).basis.n_terms - n1;
            b = zeros(n_cells*n_terms_local,CELLS(1).n_vars);
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                b(row_start:row_end,:) = CELLS(i).RHS;
            end
        end

        function b = construct_full_RHS_fully_coupled(n1,CELLS,var_idx)
            n_cells = numel(CELLS);
            n_terms_local = CELLS(1).basis.n_terms - n1;
            n_vars  = CELLS(1).n_vars;
            b = zeros(n_cells*n_terms_local*n_vars,1);
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local*n_vars + 1;
                row_end   =   i  *n_terms_local*n_vars;
                b(row_start:row_end) = CELLS(i).RHS(:);

                if (~isempty(var_idx))
                    % for n = 1:numel(CELLS(i).bc_face_id)
                    %     bc_fquad = CELLS(i).fquad( CELLS(i).bc_face_id(n) );
                    %     bc_fvec  = CELLS(i).fvec(  CELLS(i).bc_face_id(n) );
                    %     vrswt = var_rec_slip_wall_t(n1,CELLS(i),bc_fquad,bc_fvec,var_idx);
                    % 
                    %     % get the RHS addition
                    %     avgs = CELLS(i).avgs(var_idx);
                    %     tmp_RHS = vrswt.D * avgs(:);
                    %     for j = 1:numel(var_idx)
                    %         row_start = (i-1)*n_vars*n_terms_local ...
                    %                     + (var_idx(j)-1)*n_terms_local + 1;
                    %         row_end   = row_start - 1 + n_terms_local;
                    %         sub_row_start = (j-1)*n_terms_local + 1;
                    %         sub_row_end   = sub_row_start - 1 + n_terms_local;
                    %         b(row_start:row_end) = b(row_start:row_end) ...
                    %                        + tmp_RHS(sub_row_start:sub_row_end);
                    %     end
                    % end
                end
            end
        end

        function CELLS = setup_all(n1,CELLS,bc_funs)
            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            for i = 1:n_cells
                n_nbor = CELLS(i).n_nbors;
                nbor_id = zeros(n_nbor,1);
                for n = 1:n_nbor
                    nbor_id(n) = zero_mean_basis5.local_to_global( ...
                                             CELLS(i).nbor_idx(1:n_dim,n), sz );
                end
                CELLS(i) = CELLS(i).setup_rec(n1,CELLS(nbor_id),bc_funs);
            end
        end

        function CELLS = perform_iterative_reconstruction_jacobi(n1,CELLS,n_iter)
            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_vars        = CELLS(1).n_vars;
            n_cells       = numel(CELLS);

            
            for iter = 1:n_iter
                CELLS_old = CELLS;
                for i = 1:n_cells
                    n_terms = CELLS_old(i).basis.n_terms;
                    update = zeros(n_terms-n1,n_vars);
                    n_nbors = CELLS_old(i).n_nbors;
                    for n = 1:n_nbors
                        nbor_idx = CELLS_old(i).nbor_idx(1:n_dim,n);
                        j = zero_mean_basis5.local_to_global(nbor_idx,sz);
                        uj = CELLS_old(j).coefs(n1+1:n_terms,:);
                        % invA  = CELLS_old(i).self_LHS_inv;
                        AinvB = CELLS_old(i).AinvB(:,:,n);
                        % bi    = CELLS_old(i).RHS;
                        for v = 1:n_vars
                            update(:,v) = update(:,v) - AinvB * uj(:,v);
                        end
                    end
                    invA  = CELLS_old(i).self_LHS_inv;
                    bi    = CELLS_old(i).RHS;
                    for v = 1:n_vars
                        update(:,v) = update(:,v) + invA * bi(:,v);
                    end
                    CELLS(i).coefs(n1+1:n_terms,:) = update;
                end
            end
        end

        function CELLS = perform_iterative_reconstruction_SOR(n1,CELLS,omega,n_iter)
            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_vars        = CELLS(1).n_vars;
            n_cells       = numel(CELLS);
            for iter = 1:n_iter
                for i = 1:n_cells
                    n_terms = CELLS(i).basis.n_terms;
                    update = zeros(n_terms-n1,n_vars);
                    n_nbors = CELLS(i).n_nbors;
                    for n = 1:n_nbors
                        nbor_idx = CELLS(i).nbor_idx(1:n_dim,n);
                        j = zero_mean_basis5.local_to_global(nbor_idx,sz);
                        uj = CELLS(j).coefs(n1+1:n_terms,:);
                        AinvB = CELLS(i).AinvB(:,:,n);
                        for v = 1:n_vars
                            update(:,v) = update(:,v) + AinvB * uj(:,v);
                        end
                    end
                    invA  = CELLS(i).self_LHS_inv;
                    bi    = CELLS(i).RHS;
                    for v = 1:n_vars
                        update(:,v) = update(:,v) + invA * bi(:,v);
                    end
                    CELLS(i).coefs(n1+1:n_terms,:) = ...
                                  (1-omega) * CELLS(i).coefs(n1+1:n_terms,:) ...
                                   + omega  * update;
                end
            end
        end

        function [CELLS,LHS,RHS,coefs] = perform_reconstruction(n1,CELLS)
            n_terms_local = CELLS(1).basis.n_terms - n1;
            n_vars        = CELLS(1).n_vars;
            n_cells       = numel(CELLS);
            LHS = var_rec_t5.construct_full_LHS_matrix(n1,CELLS);
            RHS = var_rec_t5.construct_full_RHS(n1,CELLS);
            coefs = zeros(n_terms_local*n_cells,CELLS(1).n_vars);
            for v = 1:n_vars
                coefs(:,v) = LHS\RHS(:,v);
            end
            
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                n_terms   = CELLS(i).basis.n_terms;
                CELLS(i).coefs(n1+1:n_terms,:) = coefs(row_start:row_end,:);
            end
        end

        function [CELLS,LHS,RHS,coefs] = perform_reconstruction_fully_coupled(n1,CELLS,var_idx)
            n_terms_local = CELLS(1).basis.n_terms - n1;
            n_vars        = CELLS(1).n_vars;
            n_cells       = numel(CELLS);

            RHS = var_rec_t5.construct_full_RHS_fully_coupled(           n1, ...
                                                                      CELLS, ...
                                                                      var_idx );
            
            LHS = var_rec_t5.construct_full_LHS_matrix_fully_coupled(    n1, ...
                                                                      CELLS, ...
                                                                      var_idx );
            
            coefs = zeros(n_terms_local*n_cells,CELLS(1).n_vars);
            coefs_calc = LHS\RHS;
            for i = 1:n_cells
                n_terms   = CELLS(i).basis.n_terms;
                for j = 1:n_vars
                    row_start = (i-1)*n_vars*n_terms_local ...
                              + (j-1)*n_terms_local + 1;
                    row_end   =  row_start + n_terms_local - 1;
                    CELLS(i).coefs(n1+1:n_terms,j) = ...
                                         coefs_calc(row_start:row_end);
                end
                % for debug & comparison
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                coefs(row_start:row_end,:) = CELLS(i).coefs(n1+1:n_terms,:);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function CELLS = set_up_cell_var_recs(GRID,blk,idx_low,idx_high,degree,funs,n1,use_bcs,dist_opt)
            nvars = numel(funs);
            IDX = generate_index_array(idx_low,idx_high);
            CELLS = arrayfun( @(idx)var_rec_t5.set_up_cell(   GRID, ...
                                                               blk, ...
                                                               idx, ...
                                                            degree, ...
                                                             nvars, ...
                                                          dist_opt),   IDX );
            CELLS = arrayfun(@(CELL)CELL.set_cell_avg(funs),CELLS);
            if (use_bcs)
                CELLS = var_rec_t5.setup_all(n1,CELLS,funs);
            else
                CELLS = var_rec_t5.setup_all(n1,CELLS,[]);
            end
        end

        function CELL = set_up_cell(GRID,blk,idx,degree,n_vars,dist_opt)
            [n_dim,nbor_idx] = var_rec_t5.get_dim_and_nbor_info( GRID, ...
                                                                 blk,  ...
                                                                 [idx{:}] );

            n_nbor = size(nbor_idx,2) - 1;
            tmp_idx = num2cell( nbor_idx(:,1) );
            volume_quad = GRID.gblock(blk).grid_vars.quad( tmp_idx{:} );
            [face_quads,face_normals,nbor_face_id,bc_face_id] = var_rec_t5.get_face_info( GRID, ...
                                                                    blk, ...
                                                                    idx, ...
                                                                    n_dim );
            nodes = var_rec_t5.get_nodes(GRID,blk,nbor_idx(:,1),n_dim);
            h_ref = var_rec_t5.get_max_cell_extents(nodes,n_dim);
            basis_funs = zero_mean_basis5(n_dim,degree,h_ref,volume_quad);
            CELL = var_rec_t5( nbor_idx, nbor_face_id(1:n_nbor), bc_face_id, ...
                               basis_funs, volume_quad, face_quads, face_normals, n_vars, dist_opt );
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
            [nbor_idxs,cnt] = get_cell_face_nbor_idx( ...
                                                               idx(1:n_dim), ...
                                                                    bnd_min, ...
                                                                    bnd_max, ...
                                                                    max_dim );
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
                offset = zero_mean_basis5.global_to_local(n,nSub) - 1;
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
                if ( var_rec_t5.in_bounds(           ...
                                        nbor_idx(1:n_dim), ...
                                            ones(1,n_dim), ...
                         GRID.gblock(blk).Ncells(1:n_dim) ) )
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

        function val = in_bounds(idx,bnd_min,bnd_max)
            val = ( all(idx(:)>=bnd_min(:)) && all(idx(:)<=bnd_max(:)) || ...
                    all(idx(:)<=bnd_min(:)) && all(idx(:)>=bnd_max(:)) );
        end
    end
end