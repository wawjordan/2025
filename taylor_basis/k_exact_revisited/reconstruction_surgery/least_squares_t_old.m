classdef (Abstract) least_squares_t_old
    properties
        n_rhs              (1,1) {mustBeInteger} = 0
        n_mean_constraints (1,1) {mustBeInteger} = 0
        n_bc_constraints   (1,1) {mustBeInteger} = 0
        n_constraints      (1,1) {mustBeInteger} = 0
        LHS_m              (1,1) {mustBeInteger} = 0
        LHS_n              (1,1) {mustBeInteger} = 0
        LHS_plus           (:,:) {mustBeNumeric} = 0
        coln_scale         (:,1) {mustBeNumeric} = 0
        p                  (1,1) pascal_triangle_old
        LS_factor          (1,1) {mustBeNumeric} = 0
    end
    methods (Abstract)
        this = setup( this, n_rhs, exponents, moments, stencil, grid )
        poly = solve(this,soln_in_stencil,poly,moments,stencil,grid)
    end
    methods
        function [LHS,row] = fill_stencil_rows( this, LHS, row, col_start, exponents, moments, stencil, grid )
            % Extract sizing information
            stencil_size = size( stencil, 2 );
            n_terms = size(exponents,2);
            xi = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));
            % Fill stencil Rows in LHS
            for stencil_iter = 2:stencil_size
                % Define LHS fill indicies
                row = row + 1;
                xj = grid.gblock( stencil(1,stencil_iter) ).grid_vars.cell_c(:,stencil(2,stencil_iter),stencil(3,stencil_iter),stencil(4,stencil_iter));
                % Compute geometric weight
                LS_weight = this.get_LS_weight( xi, xj, this.LS_factor );
                for col_iter = 1:n_terms
                    % Define LHS fill indices
                    col = col_start + col_iter - 1;
                    % Compute shifted grid moment
                    moment_hat = this.compute_shifted_moment_mod_old( exponents(:,col_iter), xi, xj, exponents, moments(:,stencil_iter), this.p );
                    % moment_hat_check = this.compute_shifted_moment( exponents(:,col_iter), xi, xj, exponents, moments(:,stencil_iter), this.p );
                    % if any(abs(moment_hat - moment_hat_check)>1e-12)
                    %     [moment_hat,moment_hat_check];
                    % end

                    % Fill LHS
                    LHS(row,col) = LS_weight * moment_hat;
                end
            end
        end
    end
    methods (Static)
        function [LHS,row] = fill_mean_constraint( LHS, row, col_start, exponents, moments )
            n_terms = size(exponents,2);
            col_end = col_start + n_terms - 1;
            LHS( row, col_start:col_end ) = moments(:,1);
            row = row + 1;
        end
        function LS_weight = get_LS_weight(xi,xj,LS_factor)
            delta = xj - xi;
            LS_weight = 1/( sqrt( sum(delta.^2) )^LS_factor);
        end
        function [LHS,coln_scale] = perform_matrix_scaling(LHS)
            LHS_n = size(LHS,2);
            coln_scale = zeros(1,LHS_n);
            for n = 1:LHS_n
                coln_scale(n) = 1/sum(abs(LHS(:,n)));
            end
            for n = 1:LHS_n
                LHS(:,n) = LHS(:,n) * coln_scale(n);
            end
        end
        function [moment_hat,A,B,C] = compute_shifted_moment( pexp, xi, xj, exponents, moments, ~ )
            n_dim = numel(pexp);
            n_terms = numel(moments);

            current_exponent = zeros(3,1);
            xi_ref           = zeros(3,1);
            xj_ref           = zeros(3,1);
            current_exponent(1:n_dim) = pexp(1:n_dim);
            xi_ref(1:n_dim) = xi(1:n_dim);
            xj_ref(1:n_dim) = xj(1:n_dim);

            N = prod(current_exponent+1);
            A = zeros(N,1);
            B = zeros(N,1);
            C = zeros(N,1);

            moment_hat = 0;
            cnt = 0;
            for mp = 0:current_exponent(3)
                for mm = 0:current_exponent(2)
                    for mn = 0:current_exponent(1)
                        cnt = cnt + 1;
                        % Factorial Factor
                        factor1 = factorial(current_exponent(1))/...
                            ( factorial(mn)*factorial(current_exponent(1)-mn) );
                        factor2 = factorial(current_exponent(2))/...
                            ( factorial(mm)*factorial(current_exponent(2)-mm) );
                        factor3 = factorial(current_exponent(3))/...
                            ( factorial(mp)*factorial(current_exponent(3)-mp) );
                        factorial_factor = factor1*factor2*factor3;
                        % Grid Factor
                        factor_x = ( xj_ref(1) - xi_ref(1) )^mn;
                        factor_y = ( xj_ref(2) - xi_ref(2) )^mm;
                        factor_z = ( xj_ref(3) - xi_ref(3) )^mp;
                        grid_factor = factor_x*factor_y*factor_z;
                        % Shift Moment
                        shift_exp = current_exponent - [ mn; mm; mp ];
                        found_exp = false;
                        for shift_row = 1:n_terms
                            if ( all(exponents(1:n_dim,shift_row) == shift_exp(1:n_dim),1) )
                                found_exp = true;
                                break;
                            end
                        end

                        if ( shift_row > n_terms )
                            error(['Attempting to access entry in ',...
                                'stencil_moments which is out of bounds']);
                        end
                        if ( ~ found_exp )
                            error(['Unable to locate correct grid moment for ',...
                                'computation of shifted exponent.']);
                        end
                        shift_moment = moments(shift_row);
                        % Moment Hat
                        moment_hat = moment_hat + factorial_factor ...
                            * grid_factor ...
                            * shift_moment;
                        A(cnt) = factorial_factor;
                        B(cnt) = grid_factor;
                        C(cnt) = shift_moment;
                    end
                end
            end
        end
        function [moment_hat,A,B,C] = compute_shifted_moment_mod_old( pexp, xi, xj, exponents, moments, p )
            n_dim = numel(pexp);
            deltas = zeros(n_dim,p.Nmax);
            deltas(:,1) = 1.0;
            for m = 1:max(pexp)
                deltas(:,m+1) = (xj(1:n_dim) - xi(1:n_dim)).^m;
            end

            N = prod(pexp+1);
            A = zeros(N,1);
            B = zeros(N,1);
            C = zeros(N,1);

            moment_hat = 0;
            for m = 1:prod(pexp+1)
                cexp = least_squares_t_old.global_to_local_local(m,pexp+1)-1;
                cexp = cexp(:);

                comb_factor = 1.0;
                for j = 1:n_dim
                    comb_factor = comb_factor*p.binomial_coefficient(pexp(j),cexp(j));
                end
                grid_factor = 1.0;
                for j = 1:n_dim
                    grid_factor = grid_factor*deltas(j,cexp(j)+1);
                end
                shift_exp = pexp-cexp;
                for shift_row = 1:numel(moments)
                    if all( shift_exp == exponents(:,shift_row) )
                        break
                    end
                end
                shift_moment = moments(shift_row);

                moment_hat = moment_hat + comb_factor * grid_factor * shift_moment;
                A(m) = comb_factor;
                B(m) = grid_factor;
                C(m) = shift_moment;
            end
        end
        function iSub = global_to_local_local(iG,nSub)
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
    end

end