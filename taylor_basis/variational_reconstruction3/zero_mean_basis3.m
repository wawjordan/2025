classdef zero_mean_basis3

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
        function this = zero_mean_basis3(n_dim,degree,h_ref,quad)
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

        function x_bar = transform(this,x)
            % shift and scale the coodinate before evaluation
            scale_factor = 1.0./this.h_ref;
            x_bar = zero_mean_basis3.shift_and_scale_static( x(1:this.n_dim), this.x_ref, scale_factor );
        end

        function moment = compute_grid_moment(this,n,quad)
            % allocate
            tmp = zeros(1,quad.n_quad);

            % evaluate at quadrature points
            for q = 1:quad.n_quad
                xtmp   = this.transform( quad.quad_pts(1:this.n_dim,q) );
                % [tmp(q),coef] = this.terms(n).eval( xtmp );
                % tmp(q) = tmp(q) / coef;
                [tmp(q),~] = this.terms(n).eval( xtmp );
            end
            % integrate
            moment = quad.integrate(tmp);
        end

        function [B,fact_coef] = eval(this,n,point)
            if (n == 1) % cell avg -> piecewise constant
                B = 1.0;
                fact_coef = 1;
            else
                [B,fact_coef] = this.terms(n).eval( this.transform(point) );
                % B = B/fact_coef - this.moments(n);
                B = B - this.moments(n);
            end
        end

        function dB = deval(this,n,point,order)
            if (all(order== 0))
                dB = this.eval(n,point);
                return
            end
            if (n == 1) % cell avg -> piecewise constant
                dB = 0.0;
                return
            end
            
            [dB,dcoef,~] = this.terms(n).deval( this.transform(point), order );
            dB = dB * dcoef;
            grid_factor = 1 / prod( this.h_ref .^ order );
            dB = dB * grid_factor;
            % get normalization factor for this basis
            % [~,fact_coef1] = this.eval(n,point);
            % dB = (dB/fact_coef1) * fact_coef;
        end

        function D = calc_basis_derivative(this,term,deriv,point,scale)
            order = this.terms(deriv).exponents;
            derivative_term = this.deval( term, point, order );
            L = zero_mean_basis3.get_factorial_scaling_1(order,scale);
            D = L * derivative_term;
        end

        function D = calc_basis_normal_derivative(this,term,d_order,normal,point,scale)
            D = 0;
            start_term = this.get_n_terms(this.n_dim,d_order-1);
            end_term   = this.get_n_terms(this.n_dim,d_order);
            for i = start_term:end_term
                order = this.terms(i).exponents;
                derivative_term = this.deval( term, point, order );
                L = zero_mean_basis3.get_factorial_scaling_1(order,scale);
                Di = L * derivative_term;
                mcoeff = factorial(d_order)/prod( factorial(order) );
                D = D + mcoeff * sum(normal.^order) * Di;
            end
        end

        function deriv = calc_basis_derivatives(this,n1,point,scale)
            % Inputs:
            %  this      - zero_mean_basis2 object
            %  n1        - number of coefficients already solved for
            %  point     - point (vector) at which to evaluate the derivatives
            %  scale     - distance (scalar) describing local length scale

            % allocate the working array
            deriv = zeros(this.n_terms,this.n_terms-n1);

            % outer loop: basis functions -  n1+1:n_terms
            for j = n1+1:this.n_terms

                % inner loop: derivatives -  0:order(j-n1)
                % for i = 1:j-n1
                % for i = 1:this.n_terms
                for i = 1:j
                    % % differentiation order: [m,n,o]
                    % order = this.terms(i).exponents;
                    % 
                    % [derivative_term,~] = this.deval( j, point, order );
                    % 
                    % % delta_term = scale ^ sum( order );
                    % 
                    % L = zero_mean_basis2.get_factorial_scaling_1(order,scale);
                    % % L = scale ^ sum( order );
                    % deriv(i,j-n1) = L * derivative_term;
                    deriv(i,j-n1) = this.calc_basis_derivative(j,i,point,scale);
                end
            end
        end

        % function deriv = calc_basis_derivatives(this,n1,point,scale)
        %     % Inputs:
        %     %  this      - zero_mean_basis2 object
        %     %  n1        - number of coefficients already solved for
        %     %  point     - point (vector) at which to evaluate the derivatives
        %     %  scale     - distance (scalar) describing local length scale
        % 
        %     % allocate the working array
        %     deriv = zeros(this.n_terms,this.n_terms-n1);
        % 
        %     % outer loop: basis functions -  n1+1:n_terms
        %     for j = n1+1:this.n_terms
        % 
        %         % inner loop: derivatives -  0:order(j-n1)
        %         % for i = 1:j-n1
        %         % for i = 1:this.n_terms
        %         for i = 1:j
        %             % differentiation order: [m,n,o]
        %             order = this.terms(i).exponents;
        % 
        %             [derivative_term,~] = this.deval( j, point, order );
        % 
        %             % delta_term = scale ^ sum( order );
        % 
        %             L = zero_mean_basis2.get_factorial_scaling_1(order,scale);
        %             % L = scale ^ sum( order );
        %             deriv(i,j-n1) = L * derivative_term;
        %         end
        %     end
        % end

        function deriv = calc_basis_derivatives_vec_scale(this,n1,point,scale_vec)
            % Inputs:
            %  this      - zero_mean_basis2 object
            %  n1        - number of coefficients already solved for
            %  point     - point (vector) at which to evaluate the derivatives
            %  scale_vec - distance (vector) describing local length scale

            % allocate the working array
            deriv = zeros(this.n_terms-n1,this.n_terms-n1);

            % outer loop: basis functions -  n1+1:n_terms
            for j = n1+1:this.n_terms

                % inner loop: derivatives -  0:order(j-n1)
                for i = 1:j-n1

                    % differentiation order: [m,n,o]
                    order = this.terms(i).exponents;

                    [derivative_term,~] = this.deval( j, point, order );

                    % (dx^m dy^n dz^o ...)
                    % [delta_term,~] = this.terms(i).eval( scale_vec );
                    
                    L = zero_mean_basis3.get_factorial_scaling_2(order,scale_vec);
                    deriv(i,j-n1) = L * derivative_term;
                end
            end
        end

    end

    methods (Static)
        function L = get_factorial_scaling_1(order,scale)
            delta_term     = scale ^ sum( order );
            factorial_term = 1 / factorial( sum(order) );
            L = delta_term * factorial_term;
        end

        function L = get_factorial_scaling_2(order,scale_vec)
            n_dim = numel(order);

            delta_term = 1;
            factorial_term = 1;
            for d = 1:n_dim
                delta_term = delta_term * scale_vec(d)^order(d);
                factorial_term = factorial_term * nchoosek( sum( order(d:n_dim) ), order(d) );
            end
            den = factorial( sum(order) );
            L = delta_term * factorial_term / den;
        end

        function x_bar = shift_and_scale_static(x,shift,scale)
            x_bar = (x(:)-shift(:)).*scale(:);
        end

        function n_terms = get_n_terms( dimen, total_degree )
            n_terms = nchoosek( dimen + total_degree, total_degree );
        end

        function exp = get_exponents_sorted( dimen, total_degree )
            n_terms = zero_mean_basis3.get_n_terms(dimen,total_degree);
            exp     = zeros(dimen,n_terms);
            p_idx     = zeros(n_terms,1);
            full_degree = (total_degree+1)^dimen;
            cnt = 0;
            nSub = (total_degree+1)*ones(dimen,1);
            n10s = 10^floor(log10(full_degree));
            for n = 0:full_degree
                tmp_exp = zero_mean_basis3.global_to_local(n+1,nSub) - 1;
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