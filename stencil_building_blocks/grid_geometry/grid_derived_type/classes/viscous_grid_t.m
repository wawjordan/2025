classdef viscous_grid_t
    properties
        xi_sub_volume   (:,:,:)   {mustBeNumeric}
        eta_sub_volume  (:,:,:)   {mustBeNumeric}
        zeta_sub_volume (:,:,:)   {mustBeNumeric}
        xi_cell_area    (:,:,:,:) {mustBeNumeric}
        eta_cell_area   (:,:,:,:) {mustBeNumeric}
        zeta_cell_area  (:,:,:,:) {mustBeNumeric}
        xi_face_area    (:,:,:,:) {mustBeNumeric}
        eta_face_area   (:,:,:,:) {mustBeNumeric}
        zeta_face_area  (:,:,:,:) {mustBeNumeric}
    end
    methods
        function this = allocate_viscous_grid( this, imax, jmax, kmax, i_cells, j_cells, k_cells )
            this.xi_sub_volume   = zeros( imax   , j_cells, k_cells);
            this.eta_sub_volume  = zeros( i_cells, jmax   , k_cells);
            this.zeta_sub_volume = zeros( i_cells, j_cells, kmax   );
            this.xi_cell_area = zeros( 3, length(0:i_cells+1), j_cells            , k_cells             );
            this.xi_cell_area = zeros( 3, i_cells            , length(0:j_cells+1), k_cells             );
            this.xi_cell_area = zeros( 3, i_cells            , j_cells            , length(0:k_cells+1) );
            
            this.xi_face_area   = zeros( 3, imax, jmax, kmax );
            this.eta_face_area  = zeros( 3, imax, jmax, kmax );
            this.zeta_face_area = zeros( 3, imax, jmax, kmax );
        end
    end
end