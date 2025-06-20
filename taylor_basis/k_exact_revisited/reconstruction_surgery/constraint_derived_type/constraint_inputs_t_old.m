classdef constraint_inputs_t_old
    properties
        % Indices describing the Storage of Constrained Variables
        VAR_IDX (:,1) {mustBeInteger}

        % Evaluation & Reference Locations for the Polynomial Expansion
        x_eval  (3,1) {mustBeNumeric} = 0
        x_ref   (3,1) {mustBeNumeric} = 0

        % Unit Normal Vector
        normal  (3,1) {mustBeNumeric} = 0

        % Coefficients describing the Coupling of Different Variables
        coupling_coeff (:,1) {mustBeNumeric}
        order          (:,:) {mustBeInteger}
    end
    properties (Constant)
        NULL_CONSTRAINT_ID       =    0
        DIRICHLET_CONSTRAINT_ID  = 1000
        COUPLED_CONSTRAINT_ID    = 3000
    end
end