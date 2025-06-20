classdef constraint_t_old
    properties
        constraint_inputs    (1,1) constraint_inputs_t_old
        constraint_eval      (1,1) {mustBeNumeric} = 0
        constraint_row       (1,:) {mustBeNumeric}
        constraint_generator (1,1) constraint_generator_t_old = null_constraint_t_old(0,0,3)
    end

    methods
        function this = constraint_t_old( poly_degree, n_vars, dimen, constraint_type, constraint_eval, constraint_inputs )
            if nargin > 0
                switch constraint_type
                    case constraint_inputs_t_old.NULL_CONSTRAINT_ID
                        this.constraint_generator = null_constraint_t_old(poly_degree,n_vars,dimen);
                    case constraint_inputs_t_old.DIRICHLET_CONSTRAINT_ID
                        this.constraint_generator = dirichlet_constraint_t_old(poly_degree,n_vars,dimen);
                    case constraint_inputs_t_old.COUPLED_CONSTRAINT_ID
                        this.constraint_generator = coupled_constraint_t_old(poly_degree,n_vars,dimen);
                    otherwise
                        error('Not a valid constraint type')
                end
                this.constraint_row = zeros(1,n_vars*this.constraint_generator.nt);
                this.constraint_eval = constraint_eval;
                this.constraint_inputs = constraint_inputs;
            end
        end
        function this = generate_constraint(this,poly)
            this.constraint_row = this.constraint_generator.generate_constraint( poly, this.constraint_inputs );
        end     
    end
end