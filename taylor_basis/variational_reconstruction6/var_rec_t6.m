classdef var_rec_t6

    properties
        n_int        (1,1)   {mustBeInteger} = 1 % number of interior face nbors
        n_bnd        (1,1)   {mustBeInteger} = 0 % number of boundary face nbors
        idx_slf      (1,1)   {mustBeInteger}     % self index (global)
        idx_int      (:,1)   {mustBeInteger}     % nbor indices (global)
        f_idx_int    (:,1)   {mustBeInteger}     % corresponding face ids
        f_idx_bnd    (:,1)   {mustBeInteger}     % face ids for bcs
        basis        (1,1)   zero_mean_basis6    % local basis functions
        ord0         (1,1)   {mustBeInteger} = 0 % initial order (n in PnPm)
        n_vars       (1,1)   {mustBeInteger}     % number of eqns
        coefs        (:,:)   {mustBeNumeric}     % polynomial coeffs
        self_LHS_inv (:,:)   {mustBeNumeric}     % pseudo-inverse of self LHS
        AinvB        (:,:,:) {mustBeNumeric}     % product with self_LHS_inv
        AinvC        (:,:,:) {mustBeNumeric}     % product with self_LHS_inv
        AinvD        (:,:)   {mustBeNumeric}     % product with self_LHS_inv
        RHS          (:,:)   {mustBeNumeric}     % RHS
        get_nbor_dist
    end

    methods
        function this = var_rec_t6( idx_slf,      ...
                                    idx_int,      ...
                                    f_idx_int,    ...
                                    f_idx_bnd,    ...
                                    n_vars,       ...
                                    basis,        ...
                                    dist_opt )
            if nargin<1
                return
            end
            this.idx_slf   = idx_slf;
            this.idx_int   = idx_int;
            this.f_idx_int = f_idx_int;
            this.f_idx_bnd = f_idx_bnd;
            this.n_int     = numel(idx_int);
            this.n_vars    = n_vars;            
            this.basis     = basis;
            this.coefs     = zeros( this.basis.n_terms, n_vars );

            this.self_LHS_inv = [];
            this.AinvB        = [];
            this.AinvC        = [];
            this.AinvD        = [];
            this.RHS          = [];
            switch dist_opt
                case 'scalar_dist'
                    this.get_nbor_dist = @this.get_nbor_distance_scalar;
                case 'vector_dist'
                    this.get_nbor_dist = @this.get_nbor_distance_vec1;
                otherwise
                    error("Options are 'scalar_dist' or 'vector_dist'");
            end
        end

        function varargout = eval_reconstruction( this, x, n_terms )
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

        function this = set_cell_avg(this,fun,quad)
            n_dim   = this.basis.n_dim;
            n_quad  = quad.n_quad;
            volume  = quad.integrate(ones(1,n_quad));
            tmp_val = zeros(1,n_quad);
            for v = 1:numel(fun)
                for n = 1:n_quad
                   tmp_x      = num2cell(quad.quad_pts(1:n_dim,n));
                   tmp_val(n) = fun{v}(tmp_x{:});
                end
                this.coefs(1,v) = quad.integrate(tmp_val)/volume;
            end
        end

        % function this = setup_rec(this,n1,nbors,bc_funs)
        function this = setup_rec(this,n1,nbors,f_quad)
            this.ord0 = n1;
            self_LHS = this.get_self_LHS(nbors,f_quad);
            % if (~isempty(bc_funs))
            % end
            this.self_LHS_inv = pinv(self_LHS);
            n_terms_local = this.basis.n_terms - this.ord0;
            this.AinvB = zeros( n_terms_local, n_terms_local, this.n_int );
            nbor_LHS = this.get_nbor_LHS(nbors,f_quad);
            nbor_RHS = this.get_nbor_RHS_matrix(nbors,f_quad);
            self_RHS = this.get_self_RHS_matrix(nbors,f_quad);
            for n = 1:this.n_int
                this.AinvB(:,:,n) = this.self_LHS_inv * nbor_LHS(:,:,n);
                this.AinvC(:,:,n) = this.self_LHS_inv * nbor_RHS(:,:,n);
            end
            this.AinvD = this.self_LHS_inv * self_RHS;

            this.RHS = zeros( n_terms_local, this.n_vars );
        end

        function d = get_nbor_distance_scalar(this,nbor)
            dij = this.basis.x_ref - nbor.basis.x_ref;
            d   = sqrt( sum(dij.^2));
        end

        function d = get_nbor_distance_vec1(this,nbor)
            % d = (this.basis.h_ref + nbor.basis.h_ref)/2;
            d = abs(this.basis.x_ref - nbor.basis.x_ref);
        end

        function A = get_self_LHS(this,nbors,f_quad)
            max_n_quad = 0;
            for k = 1:this.n_int
                max_n_quad = max( max_n_quad, f_quad( this.f_idx_int(k) ).n_quad );
            end

            n_terms_local = this.basis.n_terms-this.ord0;

            d_basis = zeros( max_n_quad, this.basis.n_terms, n_terms_local );
            A = zeros( n_terms_local, n_terms_local );
            
            for k = 1:this.n_int
                dij_mag = this.get_nbor_distance_scalar(nbors(k));
                dij     = this.get_nbor_dist(nbors(k));
                n_quad  = f_quad( this.f_idx_int(k) ).n_quad;

                for q = 1:n_quad
                    point = f_quad( this.f_idx_int(k) ).quad_pts(:,q);
                    d_basis(q,:,:) = this.basis.calc_basis_derivatives( this.ord0, point, dij );
                end
                A1 = zeros( n_terms_local, n_terms_local );
                for l = 1:n_terms_local
                    for m = 1:n_terms_local
                        sum_tmp = sum( d_basis(1:n_quad,:,l) .* d_basis(1:n_quad,:,m), 2 );
                        A1(m,l) = f_quad( this.f_idx_int(k) ).integrate( sum_tmp.' );
                    end
                end
                A = A + A1/dij_mag;
            end
        end

        function D = get_self_RHS_matrix(this,nbors,f_quad)
            max_n_quad = 0;
            for k = 1:this.n_int
                max_n_quad = max( max_n_quad, f_quad( this.f_idx_int(k) ).n_quad );
            end

            n_terms_local = this.basis.n_terms - this.ord0;

            d_basis1 = zeros( max_n_quad, this.basis.n_terms, n_terms_local);
            d_basis2 = zeros( max_n_quad, this.basis.n_terms, this.ord0 );
            D = zeros( n_terms_local, this.ord0 );

            for k = 1:this.n_int
                dij_mag = this.get_nbor_distance_scalar(nbors(k));
                dij     = this.get_nbor_dist(nbors(k));
                n_quad  = f_quad( this.f_idx_int(k) ).n_quad;

                for q = 1:n_quad
                    point = f_quad( this.f_idx_int(k) ).quad_pts(:,q);
                    d_basis1(q,:,:) = this.basis.calc_basis_derivatives( this.ord0, point, dij );
                    d_basis2(q,:,:) = this.basis.calc_basis_derivatives_start( this.ord0, point, dij );
                end
                D1 = zeros( n_terms_local, this.ord0 );

                for m = 1:this.ord0
                    for l = 1:n_terms_local
                        sum_tmp = sum( d_basis1(1:n_quad,:,l) .* d_basis2(1:n_quad,:,m), 2 );
                        D1(l,m) = f_quad( this.f_idx_int(k) ).integrate( sum_tmp.' );
                    end
                end
                D = D + D1/dij_mag;
            end
        end

        function C = get_nbor_RHS_matrix(this,nbors,f_quad)
            max_n_quad = 0;
            for k = 1:this.n_int
                max_n_quad = max(max_n_quad,f_quad( this.f_idx_int(k) ).n_quad );
            end
            n_terms_local = this.basis.n_terms-this.ord0;
            d_basis1_i = zeros( max_n_quad, this.basis.n_terms, n_terms_local );
            d_basis2_j = zeros( max_n_quad, this.basis.n_terms, this.ord0 );

            C = zeros( n_terms_local, this.ord0, this.n_int );

            for k = 1:this.n_int
                dij     = this.get_nbor_dist( nbors(k) );
                dij_mag = this.get_nbor_distance_scalar( nbors(k) );
                n_quad  = f_quad( this.f_idx_int(k) ).n_quad;
                for q = 1:n_quad
                    point = f_quad( this.f_idx_int(k) ).quad_pts(:,q);
                    d_basis1_i(q,:,:) = this.basis.calc_basis_derivatives( this.ord0, point, dij );
                    d_basis2_j(q,:,:) = nbors(k).basis.calc_basis_derivatives_start( this.ord0, point, dij );
                end
                for r = 1:this.ord0
                    for l = 1:n_terms_local
                        % sum of product of derivatives
                        sum_tmp = sum( d_basis1_i(1:n_quad,:,l) .* d_basis2_j(1:n_quad,:,r), 2 );
                        C(l,r,k) = f_quad( this.f_idx_int(k) ).integrate( sum_tmp.' );
                    end
                end
                C(:,:,k) = C(:,:,k) / dij_mag;
            end
        end

        function B = get_nbor_LHS( this, nbors, f_quad )
            max_n_quad = 0;
            for k = 1:this.n_int
                max_n_quad = max( max_n_quad, f_quad( this.f_idx_int(k) ).n_quad );
            end

            n_terms_local = this.basis.n_terms - this.ord0;

            d_basis_i = zeros( max_n_quad, this.basis.n_terms, n_terms_local );
            d_basis_j = zeros( max_n_quad, this.basis.n_terms, n_terms_local );

            B = zeros( n_terms_local, n_terms_local, this.n_int );

            for k = 1:this.n_int
                dij     = this.get_nbor_dist( nbors(k) );
                dij_mag = this.get_nbor_distance_scalar( nbors(k) );
                n_quad  = f_quad( this.f_idx_int(k) ).n_quad;
                for q = 1:n_quad
                    point = f_quad( this.f_idx_int(k) ).quad_pts(:,q);
                    d_basis_i(q,:,:) = this.basis.calc_basis_derivatives( this.ord0, point, dij );
                    d_basis_j(q,:,:) = nbors(k).basis.calc_basis_derivatives( this.ord0, point, dij );
                end
                for r = 1:n_terms_local
                    for l = 1:n_terms_local
                        % sum of product of derivatives
                        sum_tmp = sum( d_basis_i(1:n_quad,:,l) .* d_basis_j(1:n_quad,:,r), 2 );
                        B(l,r,k) = f_quad( this.f_idx_int(k) ).integrate( sum_tmp.' );
                    end
                end
                B(:,:,k) = B(:,:,k) / dij_mag;
            end
        end

        function b = get_RHS_old(this,nbors,f_quad)
            n_terms_local = this.basis.n_terms - this.ord0;

            b = zeros( n_terms_local, this.n_vars );

            for k = 1:n_nbor
                dij     = this.get_nbor_dist( nbors(k) );
                dij_mag = this.get_nbor_distance_scalar( nbors(k) );
                L = zero_mean_basis6.get_length_scale( zeros(this.basis.n_dim,1), dij );
                avg_diff = nbors(k).coefs(1,:) - this.coefs(1,:);
                n_quad  = f_quad( this.f_idx_int(k) ).n_quad;
                for l = 1:n_terms_local
                    int_tmp = zeros(1,n_quad);
                    for q = 1:n_quad
                        point = f_quad( this.f_idx_int(k) ).quad_pts(:,q);
                        int_tmp(q) = this.basis.eval( l+this.ord0, point );
                    end
                    b(l,:) = b(l,:) + avg_diff(:).' * (1/dij_mag) * (L^2) * f_quad( this.f_idx_int(k) ).integrate(int_tmp);
                end
            end
        end

        function b = get_RHS(this,nbors,f_quad)
            n_terms_local = this.basis.n_terms - this.ord0;
            b = zeros( n_terms_local, this.n_vars );
            D = this.get_self_RHS_matrix(nbors,f_quad);
            C = this.get_nbor_RHS_matrix(nbors,f_quad);
            for v = 1:this.n_vars
                for k = 1:this.n_int
                    b(:,v) = b(:,v) + C(:,:,k)*nbors(k).coefs(1:this.ord0,v);
                end
                b(:,v) = b(:,v) - D*this.coefs(1:this.ord0,v);
            end
        end

        function b = get_RHS_update(this,nbors)
            n_terms_local = this.basis.n_terms - this.ord0;
            b = zeros( n_terms_local, this.n_vars );
            for v = 1:this.n_vars
                for k = 1:this.n_int
                    b(:,v) = b(:,v) + this.AinvC(:,:,k)*nbors(k).coefs(1:this.ord0,v);
                end
                b(:,v) = b(:,v) - this.AinvD*this.coefs(1:this.ord0,v);
            end
        end

    end
end