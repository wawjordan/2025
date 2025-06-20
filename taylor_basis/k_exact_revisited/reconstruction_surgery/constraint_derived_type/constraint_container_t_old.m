classdef constraint_container_t_old
    properties
        n_constraints (1,1) {mustBeInteger} = 0
        constraint    (:,1) constraint_t_old
    end
end