classdef zero_mean_basis

    properties
        degree     (1,1) {mustBeInteger} = 0    % total degree
        n_dim      (1,1) {mustBeInteger} = 0    % # of dimensions
        n_terms    (1,1) {mustBeInteger} = 0    % # of monomial terms
        terms      (:,1) monomial_t             % monomial class
        moments    (:,1) {mustBeNumeric}        % cell volume moments
        x_ref      (:,1) {mustBeNumeric}        % reference location
        h_ref      (:,1) {mustBeNumeric}        % reference length
    end

    methods
        function this = zero_mean_basis(n_dim,degree,h_ref,quad)
            % option for no arguments for initialization purposes
            if nargin < 1
                return
            end

            % set up monomials
            this.n_dim     = n_dim;
            this.degree    = degree;
            this.n_terms   = this.get_n_terms(this.n_dim,this.degree);
            exp = this.get_exponents_sorted( this.n_dim,this.degree );
            this.terms = repmat(monomial_t(),this.n_terms,1);
            for n = 1:this.n_terms
                this.terms(n).exponents = exp(:,n);
            end

            % calculate cell volume and reference lengths
            volume = quad.integrate(ones(1,quad.n_quad));
            this.x_ref     = quad.integrate(quad.quad_pts(1:this.n_dim,:))/volume;
            this.h_ref = h_ref;

            % calculate cell volume moments
            this.moments   = zeros(this.n_terms,1);
            for n = 1:this.n_terms
                this.moments(n) = this.compute_grid_moment( n, quad );
            end
            this.moments = this.moments / volume;
        end

        function [B,delta] = eval(this,n,point)
            [ delta, ~] = this.terms(n).eval( this.h_ref );
            if (n == 1)
                B = 1.0;
            else
                [B,coef] = this.terms(n).eval( this.transform(point) );
                B = B/coef - this.moments(n);
            end
        end

        function [dB,delta] = deval(this,n,point,order)
            if (all(order== 0))
                [dB,delta] = this.eval(n,point);
                return
            end
            [delta,~] = this.terms(n).deval( this.h_ref, order );
            [dB,coef] = this.terms(n).deval( this.transform(point), order );
            dB = dB/coef;
        end

        function deriv = calc_basis_derivatives_vec_scale(this,point,scale_vec)
            % allocate
            deriv = zeros(this.n_terms,this.n_terms);
            % outer loop: basis functions
            for j = 1:this.n_terms
                % inner loop: derivative
                for i = 1:j
                    order = this.terms(i).exponents;
                    [coef,~] = this.terms(i).deval( scale_vec, order );
                    [tmp1,~] = this.deval( j, point, order );
                    deriv(i,j) = tmp1 * coef;
                end
            end
        end

        function deriv = calc_basis_derivatives(this,point,scale)
            % allocate
            deriv = zeros(this.n_terms,this.n_terms);
            % outer loop: basis functions
            for j = 1:this.n_terms
                % inner loop: derivative
                for i = 1:j
                    order = this.terms(i).exponents;
                    coef = scale ^ sum( order );
                    % [ tmp1, coef1 ] = this.terms(j).deval( this.transform(point), order );
                    % deriv(i,j) = deriv(i,j) * coef / coef1;
                    [tmp1,~] = this.deval( j, point, order );
                    deriv(i,j) = tmp1 * coef;
                end
            end
        end

        function x_bar = transform(this,x)
            % shift and scale the coodinate before evaluation
            x_bar = zero_mean_basis.shift_and_scale_static( x(1:this.n_dim), this.x_ref, 1.0./this.h_ref );
        end

        function moment = compute_grid_moment(this,n,quad)
            % allocate
            tmp = zeros(1,quad.n_quad);

            % evaluate at quadrature points
            for q = 1:quad.n_quad
                xtmp   = this.transform( quad.quad_pts(1:this.n_dim,q) );
                [tmp(q),coef] = this.terms(n).deval( xtmp, 0 );
                tmp(q) = tmp(q) / coef;
            end
            % integrate
            % moment = quad.integrate(tmp) / prod( factorial( this.terms(n).exponents) );
            moment = quad.integrate(tmp);
        end

        function moment_hat = compute_shifted_moment( this, nbor, n, p )

            current_exp = this.terms(n).exponents;
            n_exp = max( current_exp );
            deltas = zeros(this.n_dim,n_exp+1);
            deltas(:,1) = 1.0;

            scales = zeros(this.n_dim,n_exp+1);
            scales(:,1) = 1.0;

            grid_factor1 = 1 / prod( factorial( current_exp ) );

            xdiff = this.transform( nbor.x_ref );
            h_ratio = nbor.h_ref./this.h_ref;

            for m = 1:n_exp
                deltas(:,m+1) = xdiff.*deltas(:,m);
                % scales(:,m+1) = ( hj(1:n_dim) ./ hi(1:n_dim) ).^m * factorial(m);
                scales(:,m+1) = h_ratio.*scales(:,m+1)*m;
            end

            moment_hat = 0;
            for m = 1:prod(current_exp+1)
                cexp = zero_mean_basis.global_to_local(m,current_exp+1)-1;

                % binomial coefficient
                comb_factor = 1.0;
                for j = 1:this.n_dim
                    comb_factor = comb_factor*p.binomial_coefficient(current_exp(j),cexp(j));
                end


                % (scaled) grid factor
                grid_factor = grid_factor1; % 1 / prod( factorial( current_exp ) )
                for j = 1:this.n_dim
                    grid_factor = grid_factor*deltas(j,cexp(j)+1);
                end

                shift_exp = current_exp-cexp(:);

                % scale factor (hj/hi)^(alpha-beta)
                scale_factor = 1.0;
                for j = 1:this.n_dim
                    scale_factor = scale_factor*scales(j,shift_exp(j)+1);
                end
                % shift moment
                for shift_row = 1:nbor.n_terms
                    if all( shift_exp == nbor.terms(shift_row).exponents, 1 )
                        break
                    end
                end
                shift_moment = nbor.moments(shift_row);
                moment_hat = moment_hat + comb_factor  ...
                                        * grid_factor  ...
                                        * scale_factor ...
                                        * shift_moment;
            end
            % get moment corresponding to current_exp in this cell
            for shift_row = 1:this.n_terms
                if all( current_exp == this.terms(shift_row).exponents, 1 )
                    break
                end
            end
            moment_hat = moment_hat - this.moments(shift_row);
        end

        function moment_hat = shift_moments(this,nbor,N,p)
            moment_hat = zeros(numel(N),1);
            for i = 1:numel(N)
                moment_hat(i) = this.compute_shifted_moment(nbor,N(i),p);
            end
        end
    end

    methods (Static)
        function x_bar = shift_and_scale_static(x,shift,scale)
            x_bar = (x(:)-shift(:)).*scale(:);
        end

        function n_terms = get_n_terms( dimen, total_degree )
            n_terms = nchoosek( dimen + total_degree, total_degree );
        end

        function exp = get_exponents_sorted( dimen, total_degree )
            n_terms = zero_mean_basis.get_n_terms(dimen,total_degree);
            exp     = zeros(dimen,n_terms);
            p_idx     = zeros(n_terms,1);
            full_degree = (total_degree+1)^dimen;
            cnt = 0;
            nSub = (total_degree+1)*ones(dimen,1);
            n10s = 10^floor(log10(full_degree));
            for n = 0:full_degree
                tmp_exp = zero_mean_basis.global_to_local(n+1,nSub) - 1;
                current_degree = sum(tmp_exp);
                if current_degree <= total_degree
                    cnt = cnt + 1;
                    exp(:,cnt) = tmp_exp;
                    % (order,reverse lexicographic in each dim) sort
                    p_idx(cnt) = current_degree*n10s*10^(dimen-1);
                    for j = 1:dimen
                        p_idx(cnt) = p_idx(cnt) + exp(j,cnt)*10^(j-1);
                    end

                end
            end
            [~,idx] = sort(p_idx);
            exp = exp(:,idx);
        end

        function iSub = global_to_local(iG,nSub)
            nDims = numel(nSub);
            iSub = zeros(1,nDims);
            if (nDims==1)
                iSub(1) = iG;
                return
            end
            p = prod(nSub);
            iGtmp = iG;
            for i = nDims:-1:1
                p = fix( p/nSub(i) );
                iTmp = mod(iGtmp-1,p)+1;
                iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
                iGtmp = iTmp;
            end
        end

        function iG = local_to_global(iSub,nSub)
            iSub = iSub(:);
            nSub = nSub(:);
            nDims = numel(iSub);
            p = 1;
            iG = 1;
            for i = 1:nDims
                iG = iG + ( iSub(i) - 1 )*p;
                p = p*nSub(i);
            end
        end

    end

end