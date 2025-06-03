classdef var_rec_t

    properties
        n_nbors  (1,1) {mustBeInteger} = 1
        idx      (:,1) {mustBeInteger}
        nbor_idx (:,:) {mustBeInteger}
        basis    (1,1) zero_mean_basis
        quad     (1,1) quad_t
        fquad    (:,1) quad_t
        fvec     (:,1) face_vec_t
        nvars    (1,1) {mustBeInteger}
        avgs     (:,1) {mustBeNumeric}
        coefs    (:,:) {mustBeNumeric}
        self_LHS (:,:) {mustBeNumeric}
        self_LHS_inv (:,:) {mustBeNumeric}
        nbor_LHS (:,:,:) {mustBeNumeric}
        RHS      (:,:) {mustBeNumeric}
    end

    methods
        function this = var_rec_t(nbor_idxs,T,Q,FQ,FN,n_vars)
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
            this.nvars    = n_vars;
            this.avgs     = zeros(n_vars,1);
            this.coefs    = zeros(this.basis.n_terms,n_vars);
            this.self_LHS = zeros(this.basis.n_terms,this.basis.n_terms);
            this.self_LHS_inv = zeros(this.basis.n_terms,this.basis.n_terms);
            this.nbor_LHS = zeros(this.basis.n_terms,this.basis.n_terms,this.n_nbors);
            this.RHS      = zeros(this.basis.n_terms,this.nvars);
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

        function this = setup_rec(this,nbors)
            this.self_LHS = this.get_self_LHS(nbors);
            this.self_LHS_inv = pinv(this.self_LHS);
            this.nbor_LHS = this.get_nbor_LHS(nbors);
            this.RHS      = this.get_RHS(nbors);
        end

        function A = get_self_LHS_CELLS_NEW(this,CELLS,idx)
            my_idx = num2cell(idx);
            n_dim = this.basis.n_dim;
            n_nbor  = this.n_nbors;
            n_terms = this.basis.n_terms;
            max_n_quad = 0;
            for k = 1:n_nbor
                max_n_quad = max(max_n_quad,this.fquad(k).n_quad );
            end
            d_basis = zeros(max_n_quad,n_terms,n_terms);
            A = zeros(n_terms,n_terms);
            for k = 1:n_nbor
                nbor_id = num2cell(CELLS(my_idx{:}).nbor_idx(1:n_dim,k));
                nbor    = CELLS(nbor_id{:});
                dij = this.basis.x_ref - nbor.basis.x_ref;
                dij_mag = sqrt( sum(dij.^2));

                n_quad = this.fquad(k).n_quad;

                for q = 1:n_quad
                    point = this.fquad(k).quad_pts(:,q);
                    % d_basis(q,:,:) = this.basis.calc_basis_derivatives(point,dij_mag);
                    d_basis(q,:,:) = this.basis.calc_basis_derivatives_vec_scale(point,dij);
                end
                for l = 1:n_terms
                    for m = 1:n_terms
                        sum_tmp = sum(d_basis(1:n_quad,:,l) .* d_basis(1:n_quad,:,m),2)/dij_mag;
                        A(m,l) = A(m,l) + this.fquad(k).integrate( sum_tmp.' );
                    end
                end
            end
        end

        function A = get_self_LHS_CELLS(this,CELLS,idx)
            my_idx = num2cell(idx);

            n_dim = this.basis.n_dim;

            n_nbor  = this.n_nbors;
            n_terms = this.basis.n_terms;
            % assume same number of quad pts on all faces
            n_quad  = this.fquad(1).n_quad;

            A = zeros(n_terms,n_terms);

            for k = 1:n_nbor
                nbor_id = num2cell(CELLS(my_idx{:}).nbor_idx(1:n_dim,k));
                nbor    = CELLS(nbor_id{:});
                dij = this.basis.x_ref - nbor.basis.x_ref;
                dij_mag = sqrt( sum(dij.^2));
                for l = 1:n_terms
                    for m = 1:n_terms
                        for q = 1:n_quad
                            point = this.fquad(k).quad_pts(:,q);
                            tmp = 0;
                            for p = 1:n_terms
                                order = this.basis.terms(p).exponents(:);
                                [tmp1,~] = this.basis.deval( l, point, order );
                                [tmp2,~] = this.basis.deval( m, point, order );
                                % [delta,~] = this.basis.terms(p).deval( dij, order );
                                % tmp = tmp + tmp1 * tmp2 * delta^2;
                                % tmp = tmp + tmp1 * tmp2 * dij_mag^(2*sum(order)-1);
                                tmp = tmp + tmp1 * tmp2 * dij_mag^(2*sum(order));
                            end
                            A(m,l) = A(m,l) + this.fquad(k).quad_wts(q)*tmp;
                        end
                    end
                end
                A = A/dij_mag;
            end
            
        end

        function A = get_self_LHS(this,nbors)
            n_nbor  = this.n_nbors;
            n_terms = this.basis.n_terms;
            max_n_quad = 0;
            for k = 1:n_nbor
                max_n_quad = max(max_n_quad,this.fquad(k).n_quad );
            end
            d_basis = zeros(max_n_quad,n_terms,n_terms);
            A = zeros(n_terms,n_terms);
            for k = 1:n_nbor
                dij = this.basis.x_ref - nbors(k).basis.x_ref;
                dij_mag = sqrt( sum(dij.^2));

                n_quad = this.fquad(k).n_quad;

                for q = 1:n_quad
                    point = this.fquad(k).quad_pts(:,q);
                    % d_basis(q,:,:) = this.basis.calc_basis_derivatives(point,dij_mag);
                    d_basis(q,:,:) = this.basis.calc_basis_derivatives_vec_scale(point,dij);
                end
                for l = 1:n_terms
                    for m = 1:n_terms
                        sum_tmp = sum(d_basis(1:n_quad,:,l) .* d_basis(1:n_quad,:,m),2)/dij_mag;
                        A(m,l) = A(m,l) + this.fquad(k).integrate( sum_tmp.' );
                    end
                end
            end
        end

        function B = get_nbor_LHS(this,nbors)
            n_nbor  = this.n_nbors;
            n_terms = this.basis.n_terms;
            max_n_quad = 0;
            for k = 1:n_nbor
                max_n_quad = max(max_n_quad,this.fquad(k).n_quad );
            end
            d_basis_1 = zeros(max_n_quad,n_terms,n_terms);
            d_basis_2 = zeros(max_n_quad,n_terms,n_terms);

            B = zeros(n_terms,n_terms,n_nbor);

            for k = 1:n_nbor
                dij = this.basis.x_ref - nbors(k).basis.x_ref;
                dij_mag = sqrt( sum(dij.^2));
                n_quad = this.fquad(k).n_quad;
                for q = 1:n_quad
                    point = this.fquad(k).quad_pts(:,q);
                    d_basis_1(q,:,:) = this.basis.calc_basis_derivatives_vec_scale(point,dij);
                    d_basis_2(q,:,:) = nbors(k).basis.calc_basis_derivatives_vec_scale(point,dij);
                end
                for l = 1:n_terms
                    for n = 1:n_terms
                        sum_tmp = sum(d_basis_1(1:n_quad,:,l) .* d_basis_2(1:n_quad,:,n),2)/dij_mag;
                        B(n,l,k) = B(n,l,k) + this.fquad(k).integrate( sum_tmp.' );
                    end
                end
            end
            
        end

        function b = get_RHS(this,nbors)
            n_nbor  = this.n_nbors;
            n_terms = this.basis.n_terms;
            % assume same number of quad pts on all faces
            n_quad  = this.fquad(1).n_quad;

            b = zeros(n_terms,1);

            for k = 1:n_nbor
                dij = this.basis.x_ref - nbors(k).basis.x_ref;
                dij_mag = sqrt( sum(dij.^2));
                for l = 1:n_terms
                    for q = 1:n_quad
                        point = this.fquad(k).quad_pts(:,q);
                        order = this.basis.terms(1).exponents(:,1);
                        b(l) = b(l) + this.fquad(k).quad_wts(q)*this.basis.deval( l, point, order ) / dij_mag;
                    end
                end
            end
        end
    end
    methods (Static)
        function A = construct_full_LHS_matrix(CELLS)

            n_dim   = CELLS(1).basis.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            n_terms = CELLS(1).basis.n_terms;

            % get the number of non-zeros
            nnz = 0;
            for i = 1:n_cells
                nnz = nnz + (CELLS(i).n_nbors + 1)*n_terms^2;
            end

            % lists for storing sparse matrix info
            % (to be used with the "sparse" command)
            row_idx = zeros(nnz,1);
            col_idx = zeros(nnz,1);
            vals    = zeros(nnz,1);

            nz_cnt = 0;
            for i = 1:n_cells
                row_start = (i-1)*n_terms + 1;
                col_start = row_start;
                [row_idx,col_idx,vals,nz_cnt] = var_rec_t.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,CELLS(i).self_LHS);
                n_nbors = CELLS(i).n_nbors;
                for n = 1:n_nbors
                    nbor_idx = CELLS(i).nbor_idx(1:n_dim,n);
                    j = zero_mean_basis.local_to_global(nbor_idx,sz);
                    col_start = (j-1)*n_terms + 1;
                    [row_idx,col_idx,vals,nz_cnt] = var_rec_t.add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,-CELLS(i).nbor_LHS(:,:,n));
                end
            end
            A = sparse(row_idx,col_idx,vals);
        end

        function [row_idx,col_idx,vals,nz_cnt] = add_array(row_idx,col_idx,vals,nz_cnt,row_start,col_start,array)
            [row_sz,col_sz] = size(array);
            cnt = nz_cnt;
            for j = 1:row_sz
                for i = 1:col_sz
                    cnt = cnt + 1;
                    row_idx(cnt) = row_start-1+i;
                    col_idx(cnt) = col_start-1+j;
                    vals(cnt)    = array(i,j);
                end
            end
            nz_cnt = cnt;
        end

        function b = construct_full_RHS(CELLS)
            n_cells = numel(CELLS);
            n_terms = CELLS(1).basis.n_terms;
            b = zeros(n_cells*n_terms,1);
            for i = 1:n_cells
                row_start = (i-1)*n_terms + 1;
                row_end   =   i  *n_terms;
                n_nbor = CELLS(i).n_nbors;
                nbors = repmat(var_rec_t(),[n_nbor,1]);
                for n = 1:n_nbor
                    idx = num2cell(CELLS(i).nbor_idx(:,n));
                    nbors(n) = CELLS(idx{:});
                end
                b(row_start:row_end) = CELLS(i).get_RHS(nbors);
            end
        end

        function CELLS = setup_all(CELLS)
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
                CELLS(i) = CELLS(i).setup_rec(CELLS(nbor_id));
            end
        end
    end
end