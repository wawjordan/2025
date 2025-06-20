classdef (Abstract) constraint_generator_t_old
    properties
        poly_order (1,1) {mustBeInteger} = 0 % polynomial order
        nv         (1,1) {mustBeInteger} = 0 % # of Variables being Reconstructed
        nt         (1,1) {mustBeInteger} = 0 % # of Terms in Polynomial Expansion
    end
    methods (Abstract)
        constraint_row = generate_constraint(this,poly,constraint_inputs)
    end
    methods
        function constraint_row = assemble_constraint(this,VAR_IDX,var_constraint)
            constraint_row = zeros(1,this.nv*this.nt);
            for vi = 1:size(VAR_IDX,1)
                % Define Extents of Constraint Contribution
                coln_start = (VAR_IDX(vi)-1)*this.nt + 1;
                coln_end   = (VAR_IDX(vi)-1)*this.nt + this.nt;

                % Error Checking
                if ( (coln_end > (this.nv * this.nt)) || (coln_start < 1) )
                    error( 'Error in constraint extraction. Attempting to fill constraint row out of bounds.' )
                end

                % Add Coefficient Contribution to Constraint Row
                constraint_row(coln_start:coln_end) = var_constraint(:,vi);
            end
        end
    end
end