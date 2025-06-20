classdef stencil_elem_mod_t_old
    properties
        elem_idx   (4,1) {mustBeInteger}
        rho        (:,1) {mustBeNumeric}
        vel        (:,1) {mustBeNumeric}
        p          (:,1) {mustBeNumeric}
        turb       (:,1) {mustBeNumeric}
        moments    (:,1) {mustBeNumeric}
        polynomial (1,1) polynomial_mod_t_old
    end
    methods
        function prim = prim(this)
            prim = [this.rho; this.vel; this.p; this.turb];
        end
    end
end