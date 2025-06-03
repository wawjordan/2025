classdef taylor_shape_functions

    properties
        degree     (1,1) {mustBeInteger} = 0    % total degree
        n_dim      (1,1) {mustBeInteger} = 0    % # of dimensions
        n_terms    (1,1) {mustBeInteger} = 0    % # of monomial terms
        exponents  (:,:) {mustBeInteger}        % monomial term exponents
        moments    (:,1) {mustBeNumeric}        % cell volume moments
        x_ref      (:,1) {mustBeNumeric}        % reference location
        h_ref      (:,1) {mustBeNumeric}        % reference length
        p          (1,1) pascal_triangle        % Pascal's triangle
    end

    methods
        function this = taylor_shape_functions(n_dim,degree,node_list,quad)
            if nargin < 1
                return
            end
            this.n_dim   = n_dim;
            this.degree  = degree;
            this.n_terms = this.get_n_terms(this.n_dim,this.degree);
            % this.exponents = this.get_exponents(this.n_dim,this.degree);
            this.exponents = this.get_exponents_sorted( this.n_dim,this.degree );
            this.h_ref     = this.get_max_cell_extents(this.n_dim,node_list);

            volume = quad.integrate(ones(1,quad.n_quad));
            this.x_ref     = quad.integrate(quad.quad_pts(1:this.n_dim,:))/volume;
            this.moments   = zeros(this.n_terms,1);
            for i = 1:this.n_terms
                this.moments(i) = this.compute_grid_moment( this.n_dim, ...
                                                            this.x_ref, ...
                                                            this.h_ref, ...
                                                            this.exponents(:,i), ...
                                                            quad );
            end
            
            this.moments = this.moments / volume;

            this.p = pascal_triangle(max(this.exponents,[],"all"));
        end
        function [B,delta] = eval(this,n,x)
            delta = this.monomial_eval(this.n_dim,this.h_ref,this.exponents(:,n));
            if (n==1)
                B = 1.0;
                return
            end
            xtmp = ( x(1:this.n_dim) - this.x_ref(:) ) ./ this.h_ref(:);
            B = this.monomial_eval(this.n_dim,xtmp,this.exponents(:,n)) ...
                / prod(factorial(this.exponents(:,n)))- this.moments(n);
        end

        function [dB,delta] = eval_derivative_alt(this,n,x,order)
            if all(order==0)
                % just evaluate the polynomial
                [dB,delta] = eval(this,n,x);
                return
            end
            exp_tmp = this.exponents(:,n) - order(:);
            delta = this.monomial_eval(this.n_dim,this.h_ref,max(exp_tmp,0));
            if (any(exp_tmp<0))
                % derivative is zero
                dB = 0.0;
                return
            end

            % find the moment with matching exponents to the differentiated monomial
            moment = 0;
            for i = 2:this.n_terms
                if all( this.exponents(:,i) == exp_tmp(:) )
                    moment = this.moments(i);
                end
            end
            xtmp = ( x(1:this.n_dim) - this.x_ref(:) ) ./ this.h_ref(:);
            dB = this.monomial_eval(this.n_dim,xtmp,exp_tmp)/prod(factorial(exp_tmp)) ...
                 - prod(factorial(this.exponents(:,n)))*moment;
        end

        function [dB,delta] = eval_derivative(this,n,x,order)
            if all(order==0)
                % just evaluate the polynomial
                [dB,delta] = eval(this,n,x);
                return
            end
            exp_tmp = this.exponents(:,n) - order(:);
            % fprintf('[%d %d]\n', exp_tmp)
            dB = 0.0;
            delta = this.monomial_eval(this.n_dim,this.h_ref,max(exp_tmp,0));
            if (any(exp_tmp<0))
                % derivative is zero
                return
            end
            xtmp = ( x(1:this.n_dim) - this.x_ref(:) ) ./ this.h_ref(:);
            dB = this.monomial_eval(this.n_dim,xtmp,exp_tmp)/prod(factorial(exp_tmp));
        end

        function val = reconstruct_from_coefs(this,coefs,x)
            val = coefs(1); % mean cell value
            for i = 2:this.n_terms
                % [B,delta] = this.eval(i,x);
                % val = val + coefs(i)*delta*B;
                [B,~] = this.eval(i,x);
                val = val + coefs(i)*B;
            end
        end

        function M = build_linear_reconstruction_matrix(this,nbor)
            N = (1+this.n_dim); % constant + linear terms
            M = zeros(1,N-1);
            M(1,:) = this.shift_moment(nbor, 2:N ) - this.moments(2:N);
        end

        function M = build_reconstruction_matrix(this,nbor)
            offset = (1+this.n_dim); % constant + linear terms
            n = this.n_terms - offset;
            M = zeros(n,n);
            M(1,:) = this.shift_moment(nbor, offset+1:this.n_terms );
            for j = 1:n
                for i = 1:n-1
                    order = this.exponents(:,i+1);
                    [dB,~] = this.eval_derivative(j+offset,nbor.x_ref,order);
                    M(i+1,j) = dB;
                end
            end
        end

        function M = build_reconstruction_matrix_alt(this,nbor)
            offset = (1+this.n_dim); % constant + linear terms
            n = this.n_terms - offset;
            M = zeros(n,n);
            M(1,:) = this.shift_moment(nbor, offset+1:this.n_terms );
            for j = 1:n
                for i = 1:n-1
                    order = this.exponents(:,i+1);
                    [dB,~] = this.eval_derivative(j+offset,nbor.x_ref,order);
                    M(i+1,j) = dB;
                end
            end
        end

        function RHS = build_linear_reconstruction_rhs(~,~,avg,nbor_avg)
            RHS = nbor_avg - avg;
        end

        function RHS = build_reconstruction_rhs(this,nbor,coefs,nbor_coefs)
            N = (1+this.n_dim); % constant + linear terms
            RHS = zeros(N,1);
            moments_hat = this.shift_moment(nbor,1:N);
            RHS(1) = nbor_coefs(1) - coefs(1);
            for i = 2:N
                RHS(1) = RHS(1) - coefs(i)*moments_hat(i);
            end

            for i = 2:N
                delta = prod( this.h_ref.^this.exponents(:,i) )/ ...
                        prod( nbor.h_ref.^this.exponents(:,i) );
                RHS(i) = delta * nbor_coefs(i) - coefs(i);
            end
        end

        function moment_hat = shift_moment(this,nbor,N)
            moment_hat = zeros(numel(N),1);
            for i = 1:numel(N)
                moment_hat(i) = taylor_shape_functions.compute_shifted_moment( ...
                                        this.exponents(:,N(i)), ...
                                        this.x_ref, this.h_ref, this.exponents, this.moments, ...
                                        nbor.x_ref, nbor.h_ref, nbor.exponents, nbor.moments, this.p );
            end
        end
    end

    methods (Static)
        function n_terms = get_n_terms( dimen, total_degree )
            n_terms = nchoosek( dimen + total_degree, total_degree );
        end

        function exp = get_exponents( dimen, total_degree )
            n_terms = taylor_shape_functions.get_n_terms(dimen,total_degree);
            exp     = zeros(dimen,n_terms);
            full_degree = (total_degree+1)^dimen;
            cnt = 0;
            nSub = (total_degree+1)*ones(dimen,1);
            for n = 0:full_degree
                tmp_exp = taylor_shape_functions.global_to_local(n+1,nSub) - 1;
                current_degree = sum(tmp_exp);
                if current_degree <= total_degree
                    cnt = cnt + 1;
                    exp(:,cnt) = tmp_exp;
                end
            end
        end

        function exp = get_exponents_sorted( dimen, total_degree )
            n_terms = taylor_shape_functions.get_n_terms(dimen,total_degree);
            exp     = zeros(dimen,n_terms);
            p_idx     = zeros(n_terms,1);
            full_degree = (total_degree+1)^dimen;
            cnt = 0;
            nSub = (total_degree+1)*ones(dimen,1);
            n10s = 10^floor(log10(full_degree));
            for n = 0:full_degree
                tmp_exp = taylor_shape_functions.global_to_local(n+1,nSub) - 1;
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

        function h_ref = get_max_cell_extents(n_dim,node_list)
            h_min = min(node_list(1:n_dim,:),[],2);
            h_max = max(node_list(1:n_dim,:),[],2);
            h_ref = 0.5 * (h_max(:) - h_min(:));
        end

        function val = monomial_eval(n_dim,x,exponents)
            val = 1.0;
            for d = 1:n_dim
                for i = exponents(d):-1:1
                    val = val*x(d);
                end
            end
        end

        function moment = compute_grid_moment(n_dim,x_ref,h_ref,exponents,quad)
            tmp = zeros(1,quad.n_quad);
            for n = 1:quad.n_quad
                xtmp   = ( quad.quad_pts(1:n_dim,n) - x_ref(1:n_dim) ) ./ h_ref(1:n_dim);
                tmp(n) = taylor_shape_functions.monomial_eval(n_dim,xtmp,exponents);
            end
            moment = quad.integrate(tmp) / prod( factorial(exponents) );
        end

        % function moment_hat = compute_shifted_moment( current_exp, xi, hi, expi, momentsi, ...
        %                                                            xj, hj, expj, momentsj )
        %     %% Setup
        %     n_termsj = numel(momentsj);
        %     moment_hat = 0;
        %     % Compute Shifted Grid Moment
        %     for mp = 0:current_exp(3)
        %         for mm = 0:current_exp(2)
        %             for mn = 0:current_exp(1)
        %                 % factorial factor
        %                 factor1 = factorial(current_exp(1))/( factorial(mn)*factorial(current_exp(1)-mn) );
        %                 factor2 = factorial(current_exp(2))/( factorial(mm)*factorial(current_exp(2)-mm) );
        %                 factor3 = factorial(current_exp(3))/( factorial(mp)*factorial(current_exp(3)-mp) );
        %                 factorial_factor = factor1*factor2*factor3;
        % 
        %                 % grid factor
        %                 factor_cx = ( (xj(1)-xi(1))/hi(1) )^mn / factorial( current_exp(1) );
        %                 factor_cy = ( (xj(2)-xi(2))/hi(2) )^mm / factorial( current_exp(2) );
        %                 factor_cz = ( (xj(3)-xi(3))/hi(3) )^mp / factorial( current_exp(3) );
        %                 grid_factor = factor_cx * factor_cy * factor_cz;
        % 
        %                 shift_exp = current_exp - [ mn; mm; mp ];
        % 
        %                 % scale factor (hj/hi)^(alpha-beta)
        %                 factor_x = ( hj(1)/hi(1) )^shift_exp(1) * factorial(shift_exp(1));
        %                 factor_y = ( hj(2)/hi(2) )^shift_exp(2) * factorial(shift_exp(2));
        %                 factor_z = ( hj(3)/hi(3) )^shift_exp(3) * factorial(shift_exp(3));
        %                 scale_factor = factor_x*factor_y*factor_z;
        % 
        %                 % shift moment
        %                 found_exp = false;
        %                 for shift_row = 1:n_termsj
        %                     if ( all(expj(:,shift_row) == shift_exp,1) )
        %                         found_exp = true;
        %                         break;
        %                     end
        %                 end
        %                 if ( shift_row > n_termsj )
        %                     error('Attempting to access entry in stencil_moments which is out of bounds');
        %                 end
        %                 if ( ~ found_exp )
        %                     error('Unable to locate correct grid moment for computation of shifted exponent.');
        %                 end
        %                 shift_moment = momentsj(shift_row);
        % 
        %                 % add contributions to moment hat
        %                 moment_hat = moment_hat + factorial_factor ...
        %                                         * grid_factor  ...
        %                                         * scale_factor     ...
        %                                         * shift_moment;
        %             end
        %         end
        %     end
        % 
        %     % get moment corresponding to current_exp in this cell
        %     n_termsi = numel(momentsi);
        %     found_exp = false;
        %     for shift_row = 1:n_termsi
        %         if ( all(expi(:,shift_row) == current_exp,1) )
        %             found_exp = true;
        %             break;
        %         end
        %     end
        %     if ( shift_row > n_termsi )
        %         error('Attempting to access entry in stencil_moments which is out of bounds');
        %     end
        %     if ( ~ found_exp )
        %         error('Unable to locate correct grid moment for computation of shifted exponent.');
        %     end
        %     self_moment = momentsi(shift_row);
        % 
        %     moment_hat = moment_hat - self_moment;
        % end
        function moment_hat = compute_shifted_moment( current_exp, xi, hi, expi, momentsi, xj, hj, expj, momentsj, p )
            n_dim = numel(current_exp);
            n_exp = max(current_exp);

            deltas = zeros(n_dim,n_exp+1);
            deltas(:,1) = 1.0;

            scales = zeros(n_dim,n_exp+1);
            scales(:,1) = 1.0;

            grid_factor1 = 1 / prod( factorial( current_exp ) );

            for m = 1:n_exp
                deltas(:,m+1) = ( (xj(1:n_dim) - xi(1:n_dim)) ./ hi(1:n_dim) ).^m;
                scales(:,m+1) = ( hj(1:n_dim) ./ hi(1:n_dim) ).^m * factorial(m);
            end

            moment_hat = 0;
            for m = 1:prod(current_exp+1)
                cexp = taylor_shape_functions.global_to_local(m,current_exp+1)-1;

                % binomial coefficient
                comb_factor = 1.0;
                for j = 1:n_dim
                    comb_factor = comb_factor*p.binomial_coefficient(current_exp(j),cexp(j));
                end


                % (scaled) grid factor
                grid_factor = grid_factor1; % 1 / prod( factorial( current_exp ) )
                for j = 1:n_dim
                    grid_factor = grid_factor*deltas(j,cexp(j)+1);
                end

                shift_exp = current_exp-cexp(:);

                % scale factor (hj/hi)^(alpha-beta)
                scale_factor = 1.0;
                for j = 1:n_dim
                    scale_factor = scale_factor*scales(j,shift_exp(j)+1);
                end
                % shift moment
                for shift_row = 1:numel(momentsj)
                    if all( shift_exp == expj(:,shift_row), 1 )
                        break
                    end
                end
                shift_moment = momentsj(shift_row);
                moment_hat = moment_hat + comb_factor  ...
                                        * grid_factor  ...
                                        * scale_factor ...
                                        * shift_moment;
            end
            % get moment corresponding to current_exp in this cell
            for shift_row = 1:numel(momentsi)
                if all( current_exp == expi(:,shift_row), 1 )
                    break
                end
            end
            moment_hat = moment_hat - momentsi(shift_row);
        end
    end

end