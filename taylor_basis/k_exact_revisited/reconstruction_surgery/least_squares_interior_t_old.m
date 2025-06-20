classdef least_squares_interior_t_old < least_squares_t_old
    methods
        function this = least_squares_interior_t_old()
            this.n_mean_constraints = 1;
            this.n_bc_constraints   = 0;
            this.n_constraints      = 1;
        end
        function this = setup( this, n_rhs, exponents, moments, stencil, grid )
            % extract sizing info
            stencil_size = size(stencil,2);
            n_terms      = size(exponents,2);

            this.n_rhs   = n_rhs;
            this.p       = pascal_triangle_old(max(exponents,[],"all"));
            
            % Set the number of rows in the LHS
            % NB: "stencil_size" does not include the current cell you are in
            stencil_size = stencil_size - 1;
            this.LHS_m = this.n_constraints + stencil_size - 1;

            % Set the number of columns in the LHS
            this.LHS_n = n_terms - 1;

            % Allocate storage for the pseudo-inverse
            this.LHS_plus = zeros( this.LHS_n, this.LHS_m );

            % Allocate storage for column scaling
            this.coln_scale = zeros(1,this.LHS_n);

            % Check Constraints
            if ( this.n_constraints > this.LHS_n )
                warning('Order of polynomial not sufficient for exactly satisfying all specified constraints.')
            end

            % Allocate Temporary Storage
            LHS = zeros( this.LHS_m, this.LHS_n );

            % Setup
            row       = 0;
            % col_start = 1;


            xi = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));
            for s = 1:stencil_size
                % Increment LHS Fill Index
                row = row + 1;

                xj = grid.gblock( stencil(1,s+1) ).grid_vars.cell_c(:,stencil(2,s+1),stencil(3,s+1),stencil(4,s+1));

                % Compute Geometric Weight
                LS_weight = this.get_LS_weight(xi,xj,this.LS_factor);
                
                for col_iter = 2:n_terms

                    % Compute shifted grid moment
                    moment_hat = this.compute_shifted_moment_mod_old( exponents(:,col_iter), xi, xj, exponents, moments(:,s+1), this.p );
                    % [moment_hat1,A1,B1,C1] = this.compute_shifted_moment_mod_old( exponents(:,col_iter), xi, xj, exponents, moments(:,s+1), this.p );
                    % [moment_hat,A2,B2,C2] = this.compute_shifted_moment( exponents(:,col_iter), xi, xj, exponents, moments(:,s+1), this.p );
                    % if any(A1~=A2)
                    %     compare = [A1,A2];
                    % end
                    % if any(abs(B1-B2)>1e-12)
                    %     compare = [B1,B2];
                    % end
                    % if any(abs(C1-C2)>1e-12)
                    %     compare = [C1,C2];
                    % end

                    % Set Column Fill Index
                    col = col_iter - 1;
                    LHS(row,col) = LS_weight*( moment_hat - moments(col_iter,1)/moments(1,1) );
                end
            end

            % scale matrix
            [LHS,this.coln_scale] = this.perform_matrix_scaling(LHS);

            % compute pseudo-inverse
            this.LHS_plus = pinv(LHS);
        end
        function poly = solve(this,soln_in_stencil,poly,moments,stencil,grid)
            % Extract Sizing Information
            % stencil_size = size( stencil, 1 )
            stencil_size = size( soln_in_stencil, 1 );
            stencil_size = stencil_size - 1;

            % Allocate Storage for RHS
            RHS = zeros( this.LHS_m, this.n_rhs );

            % Setup
            row = 0;

            % Fill Mean Information from Stencil Cells
            elem_vars = soln_in_stencil(1,:);

            xi = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));
            for si = 1:stencil_size
                % update RHS fill index
                row = row + 1;

                xj = grid.gblock( stencil(1,si+1) ).grid_vars.cell_c(:,stencil(2,si+1),stencil(3,si+1),stencil(4,si+1));

                % Compute Geometric Weight
                LS_weight = this.get_LS_weight(xi,xj,this.LS_factor);

                % fill RHS
                fit_vars = soln_in_stencil(si+1,:);
                RHS(row,:) = LS_weight * ( fit_vars - elem_vars/moments(1,1) );
            end
            % Compute reconstruction coefficients
            for row_iter = 1:this.LHS_n
                for vi = 1:this.n_rhs
                    poly.coeff(vi,row_iter+1) = this.coln_scale(row_iter) * dot( this.LHS_plus(row_iter,:),RHS(:,vi) );
                end
            end

            % Back-substitute for mean constraints
            for vi = 1:this.n_rhs
                var_mean  = elem_vars(vi);
                summation = dot( moments(2:end,1), poly.coeff(vi,2:end) );
                factor    = moments(1,1);
                poly.coeff(vi,1) = ( var_mean - summation )/factor;
            end
        end
    end
end