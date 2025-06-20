classdef dirichlet_constraint_t_old < constraint_generator_t_old
    methods
        function this = dirichlet_constraint_t_old( poly_degree, n_vars, dimen )
            % Set Sizing Information
            this.poly_order = poly_degree;
            this.nv         = n_vars;
            % this.nt          = determine_n_terms_in_polynomial( poly_degree, dimen )
            this.nt          = nchoosek(dimen+poly_degree,poly_degree);
        end
        function constraint_row = generate_constraint( this, poly, constraint_inputs )
            % Initialization
            var_constraint = zeros(this.nt,this.nv);
            % Extract Constraint
            for vi = 1:size(constraint_inputs.VAR_IDX,1)
                var_constraint(:,vi) = poly.compute_linearization( constraint_inputs.x_eval );
            end
            % Assemble constraint row
            constraint_row = this.assemble_constraint( constraint_inputs.VAR_IDX, var_constraint );
        end
    end
end