classdef least_squares_boundary_t_old < least_squares_t_old
    properties
       face_indices    (:,:) {mustBeInteger}
       bc_indices      (:,1) {mustBeInteger}
       n_boundaries    (1,1) {mustBeInteger} = 0
       constraint_eval (:,1) {mustBeNumeric}
       % same_cell       (:,1) {mustBeNumericOrLogical}
       face_id         (:,1) {mustBeInteger}
       mean_constraint_factor (1,1) {mustBeNumeric}
       boundary_constraint_factor (1,1) {mustBeNumeric}
    end
    methods
        function this = least_squares_boundary_t_old(n_var,dim,max_degree,stencil_cells,grid, mean_constraint_factor, boundary_constraint_factor )
            [n_faces,face_indices,face_id] = cell_t_old.get_face_indices(dim,max_degree,stencil_cells);
            this.n_boundaries = n_faces;
            this.face_indices = ones(4,n_faces);
            this.bc_indices   = zeros(n_faces,1);
            this.face_id      = zeros(n_faces,1);
            n_bc_constraints = 0;
            for n = 1:n_faces
                this.face_indices(1:dim+1,n) = face_indices(1:dim+1,n);
                this.face_id(n)              = face_id(n);
                b   = face_indices(1,n);
                idx = ones(1,3);
                idx(1:dim) = face_indices(2:dim+1,n);
                switch fix( ( face_id(n) + 1 )/2 )
                    case 1
                        n_quad = grid.gblock(b).grid_vars.xi_face_quad(idx(1),idx(2),idx(3)).n_quad;
                    case 2
                        n_quad = grid.gblock(b).grid_vars.eta_face_quad(idx(1),idx(2),idx(3)).n_quad;
                    case 3
                        n_quad = grid.gblock(b).grid_vars.zeta_face_quad(idx(1),idx(2),idx(3)).n_quad;
                    otherwise
                        error('invalid face id')
                end
                n_constraints = grid.gblock(b).bounds(face_id(n)).bc.n_constraints;
                n_bc_constraints = n_bc_constraints + n_quad*n_constraints;
            end
            this.n_mean_constraints = n_var;
            this.n_bc_constraints   = n_bc_constraints;
            this.n_constraints      = this.n_mean_constraints + this.n_bc_constraints;
            this.constraint_eval    = zeros(this.n_bc_constraints,1);
            % cnt = this.n_mean_constraints;
            cnt = 0;
            for n = 1:n_faces
                b   = face_indices(1,n);
                idx = ones(1,3);
                idx(1:dim) = face_indices(2:dim+1,n);
                switch fix( ( face_id(n) + 1 )/2 )
                    case 1
                        quad = grid.gblock(b).grid_vars.xi_face_quad(idx(1),idx(2),idx(3));
                    case 2
                        quad = grid.gblock(b).grid_vars.eta_face_quad(idx(1),idx(2),idx(3));
                    case 3
                        quad = grid.gblock(b).grid_vars.zeta_face_quad(idx(1),idx(2),idx(3));
                    otherwise
                        error('invalid face id')
                end
                for q = 1:quad.n_quad
                    for ci = 1:grid.gblock(b).bounds(face_id(n)).bc.n_constraints
                        cnt = cnt + 1;
                        v = grid.gblock(b).bounds(face_id(n)).VAR_IDX{ci};
                        x = quad.quad_pts(:,q);
                        if (numel(v)>1)
                            this.constraint_eval(cnt) = grid.gblock(b).bounds(face_id(n)).constraint_eval{1}(x);
                        else
                            this.constraint_eval(cnt) = grid.gblock(b).bounds(face_id(n)).constraint_eval{v}(x);
                        end
                    end
                end
            end
            this.mean_constraint_factor = mean_constraint_factor;
            this.boundary_constraint_factor = boundary_constraint_factor;
        end
        function [LHS,row,constraint_eval] = fill_bc_constraint_rows(this,dim,row,exponents,stencil,grid,LHS)
            % Extract sizing information
            % stencil_size = size( stencil, 2 );
            % n_terms      = size(exponents,2);
            poly_degree  = sum(exponents(:,end));

            n_vars       = this.n_rhs;

            idx          = stencil(2:4,1);
            x_ref        = grid.gblock(stencil(1,1)).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
            poly_local   = polynomial_mod_t_old(n_vars,dim,poly_degree);
            poly_local.x_ref(:) = x_ref(1:dim);

            constraint_eval = this.constraint_eval;

            % Setup
            bc_row = 0;
            for n = 1:this.n_boundaries
                b   = this.face_indices(1,n);
                fi = this.face_indices(2:4,n);
                switch fix( ( this.face_id(n) + 1 )/2 )
                    case 1
                        quad   = grid.gblock(b).grid_vars.xi_face_quad(fi(1),fi(2),fi(3) );
                        normal = grid.gblock(b).grid_vars.xi_n(        fi(1),fi(2),fi(3) );
                    case 2
                        quad   = grid.gblock(b).grid_vars.eta_face_quad(fi(1),fi(2),fi(3) );
                        normal = grid.gblock(b).grid_vars.eta_n(        fi(1),fi(2),fi(3) );
                    case 3
                        quad   = grid.gblock(b).grid_vars.zeta_face_quad(fi(1),fi(2),fi(3) );
                        normal = grid.gblock(b).grid_vars.zeta_n(        fi(1),fi(2),fi(3) );
                    otherwise
                        error('invalid face id')
                end
                factor = 2*mod(this.face_id(n)+1,2)-1;
                
                bc_constraints    = grid.gblock(b).bounds(n).bc;
                n_constraints     = bc_constraints.n_constraints;
                constraint_type   = grid.gblock(b).bounds(n).constraint_types;
                VAR_IDX           = grid.gblock(b).bounds(n).VAR_IDX;
                constraint_eval_1   = zeros(n_constraints,1);
                constraint_inputs = repmat(constraint_inputs_t_old(),[n_constraints,1]);
                
                for q = 1:quad.n_quad
                    for ci = 1:n_constraints
                        constraint_inputs(ci).x_eval(:)  = quad.quad_pts(:,q);
                        constraint_inputs(ci).x_ref(:)   = x_ref;
                        constraint_inputs(ci).normal(:)  = factor*normal.v(:,q);
                        constraint_inputs(ci).VAR_IDX = VAR_IDX{ci};
                        constraint_inputs(ci).coupling_coeff = factor*normal.v(:,q);
                        bc_constraints.constraint(ci) = constraint_t_old( poly_degree, n_vars, dim, constraint_type(ci), constraint_eval_1(ci), constraint_inputs(ci) );

                        % increment row counters for LHS & RHS fill
                        row = row + 1;
                        bc_row = bc_row + 1;

                        % fill LHS
                        bc_constraints.constraint(ci) = bc_constraints.constraint(ci).generate_constraint(poly_local);
                        LHS(row,:) = bc_constraints.constraint(ci).constraint_row;

                        % store RHS constraint evaluations
                        % constraint_eval(bc_row) = bc_constraints.constraint(ci).constraint_eval;
                    end
                end
            end
            % Check Constraints
            if ( row ~= this.n_constraints )
                fprintf('row = %d\n',row);
                fprintf('N constraints = %d\n',this.n_constraints);
                error('Error filling reconstruction LHS')
            end
        end
        function this = setup( this, n_rhs, exponents, moments, stencil, grid, dim )
            % extract sizing info
            stencil_size = size(stencil,2);
            n_terms      = size(exponents,2);

            this.n_rhs   = n_rhs;
            n_vars = n_rhs;
            this.p       = pascal_triangle_old(max(exponents,[],"all"));
            
            % Set the number of rows in the LHS
            % NB: "stencil_size" does not include the current cell you are in
            stencil_size = stencil_size - 1;
            this.LHS_m = this.n_constraints + ( stencil_size * n_vars );

            % Set the number of columns in the LHS
            this.LHS_n = n_terms*n_vars;

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
            row       = 1;
            col_start = 1;

            % Fill mean constraints
            for nv = 1:n_vars
                [LHS,row] = this.fill_mean_constraint(LHS,row,col_start,exponents,moments);
                col_start = col_start + n_terms;
            end

            % inconsistency - needs to be fixed
            row = row - 1;

            % Fill boundary constraints
            [LHS,row,this.constraint_eval] = fill_bc_constraint_rows(this,dim,row,exponents,stencil,grid,LHS);

            % Fill stencil rows
            col_start = 1;

            for nv = 1:n_vars
                [LHS,row] = this.fill_stencil_rows( LHS, row, col_start, exponents, moments, stencil, grid );
                col_start = col_start + n_terms;
            end

            % Error Check
            if ( row ~= this.LHS_m )
                fprintf('row = %d\n', row );
                fprintf('LHS rows %d\n', this.LHS_m );
                error( 'Error filling reconstruction LHS.' )
            end

            % Apply constraint factors to LHS
            LHS(1:n_vars,:) = this.mean_constraint_factor * LHS(1:n_vars,:);

            if ( this.n_constraints > n_vars )
                for ci = 1:this.n_bc_constraints
                    LHS(n_vars+ci,:) = this.boundary_constraint_factor * LHS(n_vars+ci,:);
                end
            end

            % scale matrix
            [LHS,this.coln_scale] = this.perform_matrix_scaling(LHS);

            % compute pseudo-inverse
            this.LHS_plus = pinv(LHS);
        end
        function poly = solve(this,soln_in_stencil,poly,~,stencil,grid)

            % Extract Sizing Information
            % stencil_size = size( stencil, 1 )
            stencil_size = size( soln_in_stencil, 1 );
            stencil_size = stencil_size - 1;

            % Allocate Storage for RHS
            RHS = zeros( this.LHS_m, 1 );

            % Setup
            n_vars = this.n_rhs;
            row    = 0;

            % Fill Mean constraint information
            fit_vars = soln_in_stencil(1,:);

            for vi = 1:n_vars

                % Update RHS fill index
                row = row + 1;

                % fill RHS
                RHS(row) = this.mean_constraint_factor * fit_vars(vi);

            end

            % fill BC constraint information
            for ci = 1:this.n_bc_constraints

                % Update RHS fill index
                row = row + 1;

                % fill RHS
                RHS(row) = this.boundary_constraint_factor * this.constraint_eval(ci);

            end

            xi = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));

            % fill information from stencil cells
            for si = 1:stencil_size

                % Update RHS fill index
                row = row + 1;

                xj = grid.gblock( stencil(1,si+1) ).grid_vars.cell_c(:,stencil(2,si+1),stencil(3,si+1),stencil(4,si+1));

                % Compute Geometric Weight
                LS_weight = this.get_LS_weight(xi,xj,this.LS_factor);

                % compute RHS
                fit_vars = soln_in_stencil(si+1,:);

                % fill RHS
                for vi = 1:n_vars
                    % Set the row offset
                    var_offset = (vi-1)*stencil_size;
                    row_offset = row + var_offset;

                    RHS(row_offset) = LS_weight * fit_vars(vi);
                end

            end

            % solve the least-squares problem
            for vi = 1:n_vars
                for ti = 1:poly.n_terms

                    % set the row number
                    row_iter = (vi-1)*poly.n_terms + ti;

                    % Compute coefficient
                    poly.coeff(vi,ti) = this.coln_scale(row_iter) * dot( this.LHS_plus(row_iter,:), RHS );
                end
            end

        end
    end
    methods (Static)
    end

end