classdef rec_cell_t
    properties
        basis       (1,1) zero_mean_basis_t
        col_scale   (:,1) {mustBeNumeric}
        coefs       (:,:) {mustBeNumeric}
        Ainv        (:,:) {mustBeNumeric}
        LHS         (:,:) {mustBeNumeric}
        n_vars      (1,1) {mustBeInteger}
        self_idx    (1,1) {mustBeInteger}
        self_block  (1,1) {mustBeInteger}
        n_nbor      (1,1) {mustBeInteger}
        n_bnd       (1,1) {mustBeInteger}
        nbor_block  (:,1) {mustBeInteger}
        nbor_idx    (:,1) {mustBeInteger}
        nbor_degree (:,1) {mustBeInteger}
        bnd_idx     (:,1) {mustBeInteger}
    end
    methods
        function this = rec_cell_t( p, self_block, self_idx, n_nbor, ...
                n_bnd, nbor_block, nbor_idx,     ...
                nbor_degree, bnd_idx, n_vars,    ...
                quad, h_ref )
            this.basis    = zero_mean_basis_t( p, quad, h_ref );
            this.n_vars   = n_vars;
            this.self_idx = self_idx;
            this.self_block = self_block;
            this.n_nbor     = n_nbor;
            this.n_bnd      = n_bnd;
            this.nbor_block = nbor_block;
            this.nbor_idx    = nbor_idx;
            this.nbor_degree = nbor_degree;
            this.bnd_idx     = bnd_idx;
            this.coefs       = zeros( p.n_terms, n_vars );
            this.LHS         = zeros( n_nbor, p.n_terms-1 );
            this.Ainv        = zeros( p.n_terms-1, n_nbor );
            this.col_scale   = zeros( p.n_terms-1, 1 );
        end

        function val = evaluate_reconstruction( this, p, point, n_terms,...
                                                n_var, var_idx )
            val(1:n_var)           = 0.0;
            local_basis(1:n_terms) = 0.0;
            for n = 1:n_terms
                local_basis(n) = this.basis.evaluate_basis( p, n, point );
            end
            for v = 1:n_var
                val(v) = val(v) + dot( this.coefs(1:n_terms,var_idx(v)),...
                                       local_basis );
            end
        end

        function dval = evaluate_reconstruction_derivative( this, p, ...
                                                            point,   ...
                                                            n_terms, ...
                                                            n_var,   ...
                                                            var_idx, ...
                                                            order )
            dval(1:n_var)          = 0.0;
            local_basis(1:n_terms) = 0.0;
            for n = 1:n_terms
                local_basis(n) = this.  ...
                                 basis. ...
                                 evaluate_basis_derivative( p,     ...
                                                            n,     ...
                                                            point, ...
                                                            order );
            end
            for v = 1:n_var
                dval(v) = dval(v) ...
                        + dot( this.coefs(1:n_terms,var_idx(v)),...
                               local_basis );
            end
        end

        function this = set_cell_avg_func( this, quad, n_dim, n_var, var_idx, eval_fun )
            tmp_val = get_cell_avg( quad, n_dim, n_var, eval_fun );
            for v = 1:n_var
                this.coefs(1,var_idx(v)) = tmp_val(v);
            end
        end

        function this = set_cell_coefs( this, n_coef, n_var, coef_idx, var_idx, vals )
            for v = 1:n_var
                for n = 1:n_coef
                    this.coefs(coef_idx(n),var_idx(v)) = vals(n,v);
                end
            end
        end
        function this = set_cell_avg_val( this, n_var, var_idx, var_val )
            this = this.set_cell_coefs( 1, n_var, 1, var_idx, var_val );
        end
        function err = get_cell_error( this, p, quad, n_terms, norm, n_var, var_idx, eval_fun )
            max_L_norm = 10;
            tmp_val = zeros(n_var,quad.n_quad);

            for n = 1:quad.nquad
                exact = eval_fun( quad.quad_pts(1:p.n_dim,n) );
                reconstructed = this.evaluate_reconstruction( p, quad.quad_pts(:,n), n_terms, n_var, var_idx );
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
        function avg = get_cell_avg( quad, n_dim, n_var, eval_fun )
            tmp_val(1:n_var,1:quad.n_quad) = 0.0;
            for n = 1:quad.n_quad
                tmp_val(:,n) = eval_fun( quad.quad_pts(1:n_dim,n) );
            end
            avg = quad.integrate( tmp_val ) / sum( quad.quad_wts );
        end
    end
end