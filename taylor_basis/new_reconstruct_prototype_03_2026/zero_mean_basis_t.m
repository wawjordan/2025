classdef zero_mean_basis_t
    properties
        moments (:,1) {mustBeNumeric}
        x_ref   (:,1) {mustBeNumeric}
        h_ref   (:,1) {mustBeNumeric}
    end

    methods
        function this = zero_mean_basis_t( p, quad, h_ref )
            if nargin > 0
                this.h_ref   = h_ref(1:p.n_dim);
                this.x_ref   = quad.                                    ...
                               integrate( quad.quad_pts(1:p.n_dim,:) )  ...
                               / sum( quad.quad_wts);
                this.moments = this.compute_grid_moments(p,quad);
            end
        end

        function x_bar = transform(this,n_dim,x)
            x_bar = (x(1:n_dim)-this.x_ref)./this.h_ref;
        end

        function moments = compute_grid_moments(this,p,quad)
            moments = zeros(p.n_terms,1);
            xtmp    = 0*quad.quad_pts(1:p.n_dim,:);
            tmp     = 0*quad.quad_wts.';
            for q = 1:quad.n_quad
                xtmp(:,q) = this.transform( p.n_dim, quad.quad_pts(:,q) );
            end
            for n = 1:p.n_terms
                for q = 1:quad.n_quad
                    [tmp(q),~] = p.eval( n, xtmp(:,q) );
                end
                moments(n) = quad.integrate(tmp);
            end
            moments = moments / sum( quad.quad_wts );
        end

        function moments = compute_shifted_moments(this,p,nbor)
            moments = zeros(p.n_terms,1);
            delta1(1:p.n_dim,1:p.total_degree+1) = 0.0;
            delta2(1:p.n_dim,1:p.total_degree+1) = 0.0;
            delta1(:,1) = 1.0;
            delta1(:,2) = this.transform( p.n_dim, nbor.x_ref );
            delta2(:,1) = 1.0;
            delta2(:,2) = nbor.h_ref ./ this.h_ref;
            for m = 2:p.total_degree
                delta1(:,m+1) = delta1(:,m) .* delta1(:,2);
                delta2(:,m+1) = delta2(:,m) .* delta2(:,2);
            end
            for n = 1:p.n_terms
                alpha = p.exponents(:,n);
                for m = 1:n
                    beta = p.exponents(:,m);
                    if (any( ( alpha - beta) < 0 ) ); continue; end
                    shift_exp = alpha - beta;
                    comb_factor = 1;
                    grid_factor = 1.0;
                    for d = 1:p.n_dim
                        comb_factor = comb_factor ...
                                    * nchoosek( alpha(d), beta(d) );
                        grid_factor = grid_factor ...
                                     * delta1( d, beta(d)+1 ) ...
                                     * delta2( d, shift_exp(d)+1 );
                    end
                    for shift_row = 1:p.n_terms
                        if ( all( p.exponents(:,shift_row) == shift_exp ) )
                            break;
                        end
                    end
                    moments(n) = moments(n) + comb_factor ...
                                            * grid_factor ...
                                            * nbor.moments(shift_row);
                end
                moments(n) = moments(n) - this.moments(n);
            end
        end
        function moments = compute_shifted_moments_quad(this,p,quad)
            moments = zeros(p.n_terms,1);
            tmp     = 0*quad.quad_wts.';
            for n = 2:p.n_terms
                for q = 1:quad.n_quad
                    tmp(q) = this.evaluate_basis(p,n,quad.quad_pts(:,q));
                end
                moments(n) = quad.integrate(tmp);
            end
            moments = moments / sum( quad.quad_wts );
        end
        function B = evaluate_basis(this,p,term,point)
            B = 1.0;
            if (term==1); return; end
            [B,~] = p.eval( term, this.transform(p.n_dim,point) );
            B = B - this.moments(term);
        end
        function dB = evaluate_basis_derivative(this,p,term,point,order)
            order = order(:);
            if ( all(order==0) )
                dB = this.evaluate_basis(p,term,point);
                return;
            end
            dB = 0.0;
            if ( term==1); return; end
            [dB,dcoef,~] = p.deval(term, this.transform(p.n_dim,point),order);
            dB = dB * dcoef / prod( this.h_ref.^order(1:p.n_dim));
        end
        function val = evaluate_reconstruction(this,p,point,coefs,n_terms,n_var,var_idx)
            val(1:n_var) = 0.0;
            basis(1:n_terms) = 0.0;
            for n = 1:n_terms
                basis(n) = this.evaluate_basis(p,n,point);
            end
            for v = 1:n_var
                val(v) = val(v) + dot(coefs(1:n_terms,var_idx(v)),basis);
            end
        end
        function dval = evaluate_reconstruction_derivative(this,p,point,coefs,n_terms,n_var,var_idx,order)
            dval(1:n_var) = 0.0;
            basis(1:n_terms) = 0.0;
            for n = 1:n_terms
                basis(n) = this.evaluate_basis_derivative(p,n,point,order);
            end
            for v = 1:n_var
                dval(v) = dval(v) + dot(coefs(1:n_terms,var_idx(v)),basis);
            end
        end
    end
end