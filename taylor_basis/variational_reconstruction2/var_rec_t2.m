classdef var_rec_t2

    properties
        n_nbors  (1,1) {mustBeInteger} = 1
        idx      (:,1) {mustBeInteger}
        nbor_idx (:,:) {mustBeInteger}
        basis    (1,1) zero_mean_basis2
        quad     (1,1) quad_t
        fquad    (:,1) quad_t
        fvec     (:,1) face_vec_t
        n_vars    (1,1) {mustBeInteger}
        avgs     (:,1) {mustBeNumeric}
        coefs    (:,:) {mustBeNumeric}
        self_LHS (:,:) {mustBeNumeric}
        self_LHS_inv (:,:) {mustBeNumeric}
        nbor_LHS (:,:,:) {mustBeNumeric}
        RHS      (:,:) {mustBeNumeric}
    end

    methods
        function this = var_rec_t2(nbor_idxs,T,Q,FQ,FN,n_vars)
            if nargin<1
                return
            end
            this.idx      = nbor_idxs(:,1);
            this.nbor_idx = nbor_idxs(:,2:end);
            this.n_nbors  = size(nbor_idxs,2) - 1;
            this.basis    = T;
            this.quad     = Q;
            this.fquad    = [FQ{:}];
            this.fvec     = [FN{:}];
            this.n_vars   = n_vars;
            this.avgs     = zeros(n_vars,1);
            this.coefs    = zeros(this.basis.n_terms,n_vars);
            this.self_LHS = zeros(this.basis.n_terms,this.basis.n_terms);
            this.self_LHS_inv = zeros(this.basis.n_terms,this.basis.n_terms);
            this.nbor_LHS = zeros(this.basis.n_terms,this.basis.n_terms,this.n_nbors);
            this.RHS      = zeros(this.basis.n_terms,this.n_vars);
        end

        function val = eval_reconstruction(this,x,n_terms)
            val = zeros(this.n_vars,1);
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
                val(v) = sum(Bx.*this.coefs(:,v));
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
        end

        function this = setup_rec(this,n1,nbors)
            this.self_LHS = this.get_self_LHS(n1,nbors);
            % this.self_LHS_inv = pinv(this.self_LHS);
            this.nbor_LHS = this.get_nbor_LHS(n1,nbors);
            this.RHS      = this.get_RHS(n1,nbors);
        end

        function A = get_self_LHS(this,n1,nbors)
            n_nbor  = this.n_nbors;
            n_terms_local = this.basis.n_terms-n1;
            max_n_quad = 0;
            for k = 1:n_nbor
                max_n_quad = max(max_n_quad,this.fquad(k).n_quad );
            end
            d_basis = zeros(max_n_quad,n_terms_local+n1,n_terms_local);
            A = zeros(n_terms_local,n_terms_local);
            
            for k = 1:n_nbor
                dij = this.basis.x_ref - nbors(k).basis.x_ref;
                % dij = abs(dij);
                dij_mag = sqrt( sum(dij.^2));
                n_quad = this.fquad(k).n_quad;

                for q = 1:n_quad
                    point = this.fquad(k).quad_pts(:,q);
                    d_basis(q,:,:) = this.basis.calc_basis_derivatives(n1,point,dij_mag);
                    % d_basis(q,:,:) = this.basis.calc_basis_derivatives_vec_scale(n1,point,dij);
                end
                A1 = zeros(n_terms_local,n_terms_local);
                for l = 1:n_terms_local
                    for m = 1:n_terms_local
                        sum_tmp = sum(d_basis(1:n_quad,:,l) .* d_basis(1:n_quad,:,m),2);
                        A1(m,l) = this.fquad(k).integrate( sum_tmp.' );
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
                max_n_quad = max(max_n_quad,this.fquad(k).n_quad );
            end
            d_basis_i = zeros(max_n_quad,n_terms_local+n1,n_terms_local);
            d_basis_j = zeros(max_n_quad,n_terms_local+n1,n_terms_local);

            B = zeros(n_terms_local,n_terms_local,n_nbor);

            for k = 1:n_nbor
                dij = this.basis.x_ref - nbors(k).basis.x_ref;
                % dij = abs(dij);
                dij_mag = sqrt( sum(dij.^2));
                n_quad = this.fquad(k).n_quad;
                for q = 1:n_quad
                    point = this.fquad(k).quad_pts(:,q);
                    % d_basis_1(q,:,:) = this.basis.calc_basis_derivatives_vec_scale(n1,point,dij);
                    % d_basis_2(q,:,:) = nbors(k).basis.calc_basis_derivatives_vec_scale(n1,point,dij);
                    d_basis_i(q,:,:) = this.basis.calc_basis_derivatives(n1,point,dij_mag);
                    d_basis_j(q,:,:) = nbors(k).basis.calc_basis_derivatives(n1,point,dij_mag);
                end
                for r = 1:n_terms_local
                    for l = 1:n_terms_local
                        % sum of product of derivatives
                        sum_tmp = sum( d_basis_i(1:n_quad,:,l) .* d_basis_j(1:n_quad,:,r), 2 );
                        B(l,r,k) = this.fquad(k).integrate( sum_tmp.' );
                    end
                end
                B(:,:,k) = B(:,:,k) / dij_mag;
            end
            
        end

        function b = get_RHS(this,n1,nbors)
            n_nbor  = this.n_nbors;
            n_terms_local = this.basis.n_terms - n1;
            % assume same number of quad pts on all faces
            n_quad  = this.fquad(1).n_quad;

            b = zeros(n_terms_local,this.n_vars);

            for k = 1:n_nbor
                dij = this.basis.x_ref - nbors(k).basis.x_ref;
                % dij = abs(dij);
                dij_mag = sqrt( sum(dij.^2));
                L = zero_mean_basis2.get_factorial_scaling_1(zeros(this.basis.n_dim,1),dij_mag);
                avg_diff = nbors(k).avgs(:) - this.avgs(:);
                for l = 1:n_terms_local
                    int_tmp = zeros(1,n_quad);
                    for q = 1:n_quad
                        point = this.fquad(k).quad_pts(:,q);
                        int_tmp(q) = this.basis.eval( l+n1, point );
                    end
                    b(l,:) = b(l,:) + avg_diff(:) * (1/dij_mag) * (L^2) * this.fquad(k).integrate(int_tmp);
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
                [row_idx,col_idx,vals,nz_cnt] = var_rec_t2.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,0*Ai);
                n_nbors = CELLS(i).n_nbors;
                for n = 1:n_nbors
                    nbor_idx = CELLS(i).nbor_idx(1:n_dim,n);
                    j = zero_mean_basis.local_to_global(nbor_idx,sz);
                    col_start = (j-1)*n_terms_local + 1;
                    Bij = CELLS(i).nbor_LHS(:,:,n);
                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t2.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-Bij);
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
                n_nbor = CELLS(i).n_nbors;
                nbors = repmat(var_rec_t2(),[n_nbor,1]);
                for n = 1:n_nbor
                    idx = num2cell(CELLS(i).nbor_idx(:,n));
                    nbors(n) = CELLS(idx{:});
                end
                bi = CELLS(i).get_RHS(n1,nbors);
                b(row_start:row_end,:) = bi;
            end
        end

        function CELLS = setup_all(n1,CELLS)
            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            for i = 1:n_cells
                n_nbor = CELLS(i).n_nbors;
                nbor_id = zeros(n_nbor,1);
                for n = 1:n_nbor
                    nbor_id(n) = zero_mean_basis.local_to_global(CELLS(i).nbor_idx(1:n_dim,n),sz);
                end
                CELLS(i) = CELLS(i).setup_rec(n1,CELLS(nbor_id));
            end
        end

        function [CELLS,LHS,RHS,coefs] = perform_reconstruction(n1,CELLS)
            n_terms_local = CELLS(1).basis.n_terms - n1;
            n_vars        = CELLS(1).n_vars;
            n_cells       = numel(CELLS);
            LHS = var_rec_t2.construct_full_LHS_matrix(n1,CELLS);
            RHS = var_rec_t2.construct_full_RHS(n1,CELLS);
            coefs = zeros(n_terms_local*n_cells,CELLS(1).n_vars);
            for v = 1:n_vars
                coefs(:,v) = LHS\RHS(:,v);
            end
            
            for i = 1:n_cells
                row_start = (i-1)*n_terms_local + 1;
                row_end   =   i  *n_terms_local;
                n_terms   = CELLS(i).basis.n_terms;
                CELLS(i).coefs(1,:) = CELLS(i).avgs;
                CELLS(i).coefs(n1+1:n_terms,:) = coefs(row_start:row_end,:);
            end
        end
    end
end