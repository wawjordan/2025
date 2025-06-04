classdef cell_soln_info

    properties
        n_nbors  (1,1) {mustBeInteger} = 1
        idx      (:,1) {mustBeInteger}
        nbor_idx (:,:) {mustBeInteger}
        avgs     (:,1) {mustBeNumeric}
        rec_LHS  (:,1) cell
        rec_RHS  (:,1) cell
        taylor   (1,1) taylor_shape_functions
        quad     (1,1) quad_t
        fquad    (:,1) quad_t
        fvec     (:,1) face_vec_t
    end

    methods
        function this = cell_soln_info(nbor_idxs,T,Q,FQ,FN,n_vars)
            this.taylor   = T;
            this.quad     = Q;
            this.idx      = nbor_idxs(:,1);
            this.nbor_idx = nbor_idxs(:,2:end);
            this.n_nbors  = size(nbor_idxs,2) - 1;
            this.fquad    = [FQ{:}];
            this.fvec     = [FN{:}];
            this.avgs     = zeros(n_vars,1);
            degree = this.taylor.degree;
            n_dim  = this.taylor.n_dim;
            this.rec_LHS = cell(degree,1);
            this.rec_RHS = cell(degree,1);

            % n_terms = 1; % 0th degree - cell average
            for i = 1:degree
                n_terms = this.taylor.get_n_terms(n_dim,i) - 1;
                this.rec_LHS{i} = zeros(n_terms,this.n_nbors);
                this.rec_RHS{i} = zeros(n_terms,this.n_nbors);
            end
        end
        function this = set_cell_avg(this,funs)
            n_dim = this.taylor.n_dim;
            for i = 1:numel(this.avgs)
                this.avgs(i) = cell_soln_info.get_cell_avg(this.quad,funs{i},n_dim);
            end
        end
        function this = set_rec_LHS1(this,CELLS)
            N = (1+this.taylor.n_dim); % constant + linear terms
            for j = 1:this.n_nbors
                tmp_idx = num2cell(this.nbor_idx(:,j));
                this.rec_LHS{1}(:,j) = this.taylor.shift_moment( CELLS(tmp_idx{:}).taylor, 2:N );
            end
        end
        function this = set_rec_LHS(this,CELLS)
            for j = 1:this.n_nbors
                tmp_idx = num2cell(this.nbor_idx(:,j));
                this.rec_LHS{1}(j,:) = this.taylor.build_linear_reconstruction_matrix( ...
                    CELLS(tmp_idx{:}).taylor );
            end
        end
    end
    methods (Static)
        function avg = get_cell_avg(Q,fun,n_dim)
            tmp_val = zeros(1,Q.n_quad);
            for n = 1:Q.n_quad
                tmp_x      = num2cell(Q.quad_pts(1:n_dim,n));
                tmp_val(n) = fun(tmp_x{:});
            end
            volume = Q.integrate(ones(1,Q.n_quad));
            avg    = Q.integrate(tmp_val)/volume;
        end
        function dvals = get_derivative_approximation(CELLS,idx)
            my_idx = num2cell(idx);
            n_nbors = CELLS(my_idx{:}).n_nbors;
            n_dim = CELLS(my_idx{:}).taylor.n_dim;
            M = zeros(n_nbors-1,n_dim);
            R = zeros(n_nbors-1,1);

            n_vars = numel(CELLS(my_idx{:}).avgs);
            dvals = zeros(n_dim,n_vars);
            for n = 1:n_vars
                for j = 2:n_nbors
                    nbor_idx = num2cell(CELLS(my_idx{:}).nbor_idx(:,j));
                    M(j-1,:) = CELLS(my_idx{:}).taylor.build_linear_reconstruction_matrix( ...
                        CELLS(nbor_idx{:}).taylor );
                    % [CELLS(my_idx{:}).taylor.moments,CELLS(nbor_idx{:}).taylor.moments]
                    R(j-1)   = CELLS(my_idx{:}).taylor.build_linear_reconstruction_rhs( ...
                        CELLS(nbor_idx{:}).taylor, ...
                        CELLS(my_idx{:}).avgs(n), ...
                        CELLS(nbor_idx{:}).avgs(n) );
                end
                % M
                % R
                dvals(:,n) = M\R;
            end
        end
        function dvals = get_derivative_approximation_moments(CELLS,idx)
            my_idx = num2cell(idx);
            n_nbors = CELLS(my_idx{:}).n_nbors;
            n_dim = CELLS(my_idx{:}).taylor.n_dim;
            M = zeros(n_nbors-1,n_dim);
            R = zeros(n_nbors-1,1);

            n_vars = 2;
            dvals = zeros(n_dim,n_vars);
            for n = 1:n_vars
                for j = 2:n_nbors
                    nbor_idx = num2cell(CELLS(my_idx{:}).nbor_idx(:,j));
                    M(j-1,:) = CELLS(my_idx{:}).taylor.build_linear_reconstruction_matrix( ...
                        CELLS(nbor_idx{:}).taylor );
                    % [CELLS(my_idx{:}).taylor.moments,CELLS(nbor_idx{:}).taylor.moments]
                    R(j-1)   = CELLS(my_idx{:}).taylor.build_linear_reconstruction_rhs( ...
                        CELLS(nbor_idx{:}).taylor, ...
                        CELLS(my_idx{:}).taylor.moments(n+1), ...
                        CELLS(nbor_idx{:}).taylor.moments(n+1) );
                end
                % M
                % R
                dvals(:,n) = M\R;
            end
        end

        function face_deriv = calc_basis_derivatives(CELL,point)
            n_terms = CELL.taylor.n_terms;
            face_deriv = zeros(n_terms,n_terms);

            % outer loop: basis functions
            for j = 1:n_terms
                % inner loop: derivative
                for i = 1:j
                    order = CELL.taylor.exponents(:,i);
                    face_deriv(i,j) = CELL.taylor.eval_derivative(j,point,order);
                end
            end
        end

        function face_deriv = calc_basis_derivatives_scaled(CELL,point,dij)
            n_terms = CELL.taylor.n_terms;
            face_deriv = zeros(n_terms,n_terms);

            % outer loop: basis functions
            for j = 1:n_terms
                % inner loop: derivative
                for i = 1:j
                    order = CELL.taylor.exponents(:,i);
                    face_deriv(i,j) = CELL.taylor.eval_derivative(j,point,order) * dij^( sum(order) );
                end
            end
        end

        function [A,B,b] = calc_cell_matrices(CELLS,idx)
            my_idx = num2cell(idx);
            CELL = CELLS(my_idx{:});
            n_nbors = CELL.n_nbors;
            n_terms = CELL.taylor.n_terms;

            A = zeros(n_terms,n_terms);
            B = zeros(n_terms,n_terms,n_nbors);
            b = zeros(n_terms,1);

            % assume same number of quad pts on all faces
            n_quad  = CELL.fquad(1).n_quad;

            d_basis_1 = zeros(n_quad,n_terms,n_terms);
            d_basis_2 = zeros(n_quad,n_terms,n_terms);

            for k = 1:n_nbors

                % get neighbor
                nbor_idx = num2cell(CELL.nbor_idx(:,k));
                NBOR = CELLS(nbor_idx{:});

                dij = sqrt( sum( CELL.taylor.x_ref - NBOR.taylor.x_ref).^2 );
                % dij = 1;

                % collect basis function derivatives in working arrays
                for q = 1:n_quad
                    point = CELL.fquad(k).quad_pts(:,q);
                    % d_basis_1(q,:,:) = cell_soln_info.calc_basis_derivatives(CELL,point);
                    % d_basis_2(q,:,:) = cell_soln_info.calc_basis_derivatives(NBOR,point);
                    d_basis_1(q,:,:) = cell_soln_info.calc_basis_derivatives_scaled(CELL,point,dij);
                    d_basis_2(q,:,:) = cell_soln_info.calc_basis_derivatives_scaled(NBOR,point,dij);
                end

                for l = 1:n_terms
                    for m = 1:n_terms
                        A(m,l) = A(m,l) + CELL.fquad(k).integrate( sum(d_basis_1(:,:,l) .* d_basis_1(:,:,m),2).' ) / dij;
                    end

                    for n = 1:n_terms
                        B(n,l,k) = B(n,l,k) + CELL.fquad(k).integrate( sum( d_basis_1(:,:,l) .* d_basis_2(:,:,n),2).' ) / dij;
                    end
                    b(l) = b(l) + CELL.fquad(k).integrate( d_basis_1(:,1,l).' ) / dij;
                end
            end

        end

        function A = construct_full_LHS_matrix(CELLS)

            n_dim   = CELLS(1).taylor.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            n_terms = CELLS(1).taylor.n_terms;

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
                idx = CELLS(1).taylor.global_to_local(i,sz);
                [D,Bij,~] = CELLS(1).calc_cell_matrices_explicit(CELLS,idx);

                row_start = (i-1)*n_terms + 1;
                col_start = row_start;
                [row_idx,col_idx,vals,nz_cnt] = CELLS(1).add_array(row_idx,col_idx,vals,nz_cnt,...
                    row_start,col_start,D);
                n_nbors = CELLS(i).n_nbors;
                for n = 1:n_nbors
                    nbor_idx = CELLS(i).nbor_idx(1:n_dim,n);
                    j = CELLS(1).taylor.local_to_global(nbor_idx,sz);
                    col_start = (j-1)*n_terms + 1;
                    [row_idx,col_idx,vals,nz_cnt] = CELLS(1).add_array(row_idx,col_idx,vals,nz_cnt,...
                        row_start,col_start,-Bij(:,:,n));
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
                    vals(cnt)     = array(i,j);
                end
            end
            nz_cnt = cnt;
        end


        function [A,B,b] = calc_cell_matrices_explicit(CELLS,idx)
            my_idx = num2cell(idx);
            CELL = CELLS(my_idx{:});
            n_nbors = CELL.n_nbors;
            n_terms = CELL.taylor.n_terms;

            A = zeros(n_terms,n_terms);
            B = zeros(n_terms,n_terms,n_nbors);
            b = zeros(n_terms,1);

            % assume same number of quad pts on all faces
            n_quad  = CELL.fquad(1).n_quad;
            cnt = 0;
            for k = 1:n_nbors

                % get neighbor
                nbor_idx = num2cell(CELL.nbor_idx(:,k));
                NBOR = CELLS(nbor_idx{:});

                dij = sqrt( sum( CELL.taylor.x_ref - NBOR.taylor.x_ref).^2 );
                % dij = 1;
                for l = 1:n_terms
                    for m = 1:n_terms
                        for q = 1:n_quad
                            point = CELL.fquad(k).quad_pts(:,q);
                            tmp = 0;
                            for p = 1:n_terms
                                order = CELL.taylor.exponents(:,p);
                                tmp1 = CELL.taylor.eval_derivative(l,point,order);
                                tmp2 = CELL.taylor.eval_derivative(m,point,order);
                                tmp = tmp + tmp1 * tmp2 * dij^(2*sum(order)-1);
                            end
                            A(m,l) = A(m,l) + CELL.fquad(k).quad_wts(q)*tmp;
                        end
                    end

                    for n = 1:n_terms
                        for q = 1:n_quad
                            point = CELL.fquad(k).quad_pts(:,q);
                            tmp = 0;
                            for p = 1:n_terms
                                order = CELL.taylor.exponents(:,p);
                                tmp = tmp + CELL.taylor.eval_derivative(l,point,order) ...
                                    * NBOR.taylor.eval_derivative(n,point,order) ...
                                    * dij^(2*sum(order)-1);
                            end
                            B(n,l,k) = B(n,l,k) + CELL.fquad(k).quad_wts(q)*tmp;
                        end
                    end
                    for q = 1:n_quad
                        point = CELL.fquad(k).quad_pts(:,q);
                        % b(l) = b(l) + CELL.fquad(k).quad_wts(q)*CELL.taylor.eval(l,point);
                        order = CELL.taylor.exponents(:,1);
                        b(l) = b(l) + CELL.fquad(k).quad_wts(q)*CELL.taylor.eval_derivative(l,point,order) / dij;
                    end
                end
            end

        end

        function b = construct_full_RHS(CELLS)

            n_dim   = CELLS(1).taylor.n_dim;
            dim_sz  = num2cell(1:n_dim);
            sz      = size(CELLS,dim_sz{:});
            n_cells = numel(CELLS);
            n_terms = CELLS(1).taylor.n_terms;

            b = zeros(n_cells*n_terms,1);
            for i = 1:n_cells
                idx = CELLS(1).taylor.global_to_local(i,sz);
                row_start = (i-1)*n_terms + 1;
                row_end   =   i  *n_terms;
                b(row_start:row_end) = CELLS(1).calc_cell_rhs_explicit(CELLS,idx);
            end

        end

        function b = calc_cell_rhs_explicit(CELLS,idx)

            % currently hard-coded
            var = 1;

            my_idx = num2cell(idx);
            CELL = CELLS(my_idx{:});
            n_nbors = CELL.n_nbors;
            n_terms = CELL.taylor.n_terms;

            b = zeros(n_terms,1);

            % assume same number of quad pts on all faces
            n_quad  = CELL.fquad(1).n_quad;

            for k = 1:n_nbors

                % get neighbor
                nbor_idx = num2cell(CELL.nbor_idx(:,k));
                NBOR = CELLS(nbor_idx{:});
                dij = sqrt( sum( CELL.taylor.x_ref - NBOR.taylor.x_ref).^2 );
                % dij = 1;

                for l = 1:n_terms
                    for q = 1:n_quad
                        point = CELL.fquad(k).quad_pts(:,q);
                        % b(l) = b(l) + CELL.fquad(k).quad_wts(q)*CELL.taylor.eval(l,point);
                        order = CELL.taylor.exponents(:,1);
                        avg_diff = NBOR.avgs(var) - CELL.avgs(var);
                        b(l) = b(l) + avg_diff * CELL.fquad(k).quad_wts(q)*CELL.taylor.eval_derivative(l,point,order) / dij;
                    end
                end
            end

        end

    end
end
