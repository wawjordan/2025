classdef reconstruct_mod_t_old
    properties
        stencil_size  (1,1) {mustBeInteger}
        stencil       (:,:) {mustBeInteger}
        polynomial    (1,1) polynomial_mod_t_old
        moments       (:,:) {mustBeNumeric}
        least_squares (1,1) least_squares_t_old = least_squares_interior_t_old()
    end
    methods
        function this = reconstruct_mod_t_old()
        end
        function this = setup_interior( this, n_var, n_dim, degree, stencil, grid )
            % Extract sizing info
            this.stencil_size = size(stencil,2);

            % "create" reconstruction stencils
            this.stencil      = stencil;
            x_ref = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));

            % Setup reconstruction polynomial
            this.polynomial = polynomial_mod_t_old( n_var, n_dim, degree );
            this.polynomial.x_ref = x_ref(:);

            % allocate grid moments
            this.moments = zeros( this.polynomial.n_terms, this.stencil_size );

            % compute grid moments
            for i = 1:this.stencil_size
                x_ref = grid.gblock( stencil(1,i) ).grid_vars.cell_c(:,stencil(2,i),stencil(3,i),stencil(4,i));
                quad  = grid.gblock( stencil(1,i) ).grid_vars.quad(stencil(2,i),stencil(3,i),stencil(4,i));
                this.moments(:,i) = compute_grid_moments_mod_old( this.polynomial.exponents, x_ref(:), quad );
            end
            % this = setup( this, n_rhs, exponents, moments, stencil, grid )
            this.least_squares = least_squares_interior_t_old();
            this.least_squares = this.least_squares.setup( n_var, this.polynomial.exponents, this.moments, stencil, grid );
            % this.least_squares = least_squares_interior_t_old( n_var, this.polynomial.exponents, this.moments, stencil, grid, type );
        end
        function this = setup_boundary( this, n_var, n_dim, degree, stencil, grid, neighbor_max_degree, stencil_cells, mean_constraint_factor, boundary_constraint_factor )
            
            % Extract sizing info
            this.stencil_size = size(stencil,2);

            % "create" reconstruction stencils
            this.stencil      = stencil;
            x_ref = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));

            % Setup reconstruction polynomial
            this.polynomial = polynomial_mod_t_old( n_var, n_dim, degree );
            this.polynomial.x_ref = x_ref(:);

            % allocate grid moments
            this.moments = zeros( this.polynomial.n_terms, this.stencil_size );

            % compute grid moments
            for i = 1:this.stencil_size
                x_ref = grid.gblock( stencil(1,i) ).grid_vars.cell_c(:,stencil(2,i),stencil(3,i),stencil(4,i));
                quad  = grid.gblock( stencil(1,i) ).grid_vars.quad(stencil(2,i),stencil(3,i),stencil(4,i));
                this.moments(:,i) = compute_grid_moments_mod_old( this.polynomial.exponents, x_ref(:), quad );
            end
            this.least_squares = least_squares_boundary_t_old( n_var, n_dim, neighbor_max_degree, stencil_cells, grid, mean_constraint_factor, boundary_constraint_factor);
            this.least_squares  = this.least_squares.setup( n_var, this.polynomial.exponents, this.moments, stencil, grid, n_dim );
        end
        function this = setup_boundary2( this, n_var, n_dim, degree, stencil, grid, neighbor_max_degree, stencil_cells )
            
            % Extract sizing info
            this.stencil_size = size(stencil,2);

            % "create" reconstruction stencils
            this.stencil      = stencil;
            x_ref = grid.gblock( stencil(1,1) ).grid_vars.cell_c(:,stencil(2,1),stencil(3,1),stencil(4,1));

            % Setup reconstruction polynomial
            this.polynomial = polynomial_mod_t_old( n_var, n_dim, degree );
            this.polynomial.x_ref = x_ref(:);

            % allocate grid moments
            this.moments = zeros( this.polynomial.n_terms, this.stencil_size );

            % compute grid moments
            for i = 1:this.stencil_size
                x_ref = grid.gblock( stencil(1,i) ).grid_vars.cell_c(:,stencil(2,i),stencil(3,i),stencil(4,i));
                quad  = grid.gblock( stencil(1,i) ).grid_vars.quad(stencil(2,i),stencil(3,i),stencil(4,i));
                this.moments(:,i) = compute_grid_moments_mod_old( this.polynomial.exponents, x_ref(:), quad );
            end
            this.least_squares = least_squares_boundary2_t_old( n_var, n_dim, neighbor_max_degree, stencil_cells, grid);
            this.least_squares  = this.least_squares.setup( n_var, this.polynomial.exponents, this.moments, stencil, grid, n_dim );
        end
    end
end