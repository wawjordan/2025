classdef derived_grid_vars_t_old
    properties
        % centroidal coordinates for each cell of each block
        cell_c (:,:,:,:) {mustBeNumeric}

        xi_area   (:,:,:) {mustBeNumeric}
        eta_area  (:,:,:) {mustBeNumeric}
        zeta_area (:,:,:) {mustBeNumeric}
        volume    (:,:,:) {mustBeNumeric}
        % face normals (3,imax,jmax,kmax)
        % xi_n   (:,:,:,:) 
        % eta_n  (:,:,:,:) 
        % zeta_n (:,:,:,:)
        xi_n   (:,:,:) face_vec_t_old
        eta_n  (:,:,:) face_vec_t_old
        zeta_n (:,:,:) face_vec_t_old
        % volume quadrature information (icells, jcells, kcells)
        quad (:,:,:) quad_t_old
        % face quadrature information
        xi_face_quad   (:,:,:) quad_t_old
        eta_face_quad  (:,:,:) quad_t_old
        zeta_face_quad (:,:,:) quad_t_old 
    end
    methods
        function this = allocate_derived_grid( this, i_cells, j_cells, k_cells )
            this.volume = zeros(    i_cells, j_cells, k_cells );
            this.cell_c = zeros( 3, i_cells, j_cells, k_cells );
            this.xi_area   = zeros(    i_cells+1, j_cells  , k_cells   );
            this.eta_area  = zeros(    i_cells  , j_cells+1, k_cells   );
            this.zeta_area = zeros(    i_cells  , j_cells  , k_cells+1 );
            % this.xi_n      = zeros( 3, i_cells+1, j_cells  , k_cells   );
            % this.eta_n     = zeros( 3, i_cells  , j_cells+1, k_cells   );
            % this.zeta_n    = zeros( 3, i_cells  , j_cells  , k_cells+1 );
            this.xi_n   = repmat(face_vec_t_old(), i_cells+1, j_cells  , k_cells   );
            this.eta_n  = repmat(face_vec_t_old(), i_cells  , j_cells+1, k_cells   );
            this.zeta_n = repmat(face_vec_t_old(), i_cells  , j_cells  , k_cells+1 );
            
            this.quad           = repmat(quad_t_old(), i_cells  , j_cells  , k_cells   );
            this.xi_face_quad   = repmat(quad_t_old(), i_cells+1, j_cells  , k_cells   );
            this.eta_face_quad  = repmat(quad_t_old(), i_cells  , j_cells+1, k_cells   );
            this.zeta_face_quad = repmat(quad_t_old(), i_cells  , j_cells  , k_cells+1 );
        end
    end
end