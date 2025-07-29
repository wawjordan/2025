classdef var_rec_slip_wall_t
    % container class for routines handling slip wall bc
properties
        n_mtm       (1,1)   {mustBeInteger}     % number of eqns
        mtm_idx     (:,1)   {mustBeInteger}     % index in variable vector
        % EinvD       (:,:)   {mustBeNumeric}
        E           (:,:)   {mustBeNumeric}
        D           (:,:)   {mustBeNumeric}
    end

    methods
        function this = var_rec_slip_wall_t(n1,var_rec,bc_fquad,bc_fvec,var_idx)
            if nargin<1
                return
            end
            this.n_mtm = numel(var_idx);
            [this.E,this.D] = this.get_E_and_D_matrix(n1,var_rec,bc_fquad,bc_fvec);
            % this.EinvD    = zeros(this.basis.n_terms,this.basis.n_terms,this.n_nbors);
            % this.EinvD = E\D; % inv(E)*D
        end

        function [E,D] = get_E_and_D_matrix(this,n1,var_rec,bc_fquad,bc_fvec)

            % number of unsolved coefficients
            n_terms_local = var_rec.basis.n_terms-n1;

            % number of quadrature points for this boundary interface
            n_quad = bc_fquad.n_quad;

            % area of the boundary face
            area = bc_fquad.integrate(ones(1,n_quad));

            % area centroid of the boundary face
            face_x_ref = bc_fquad.integrate(bc_fquad.quad_pts(1:var_rec.basis.n_dim,:))/area;

            % distance between cell centroid and boundary face centroid
            dbf = var_rec.basis.x_ref - face_x_ref;

            % magnitude of dbf
            dbf_mag = sqrt( sum(dbf.^2) );

            % array containing various derivatives of the basis functions
            % d_basis = zeros(n_quad,var_rec.basis.degree+1,n_terms_local);
            d_basis = zeros(n_quad,n_terms_local+1,n_terms_local);
            for q = 1:n_quad
                point  = bc_fquad.quad_pts(:,q);
                normal = bc_fvec.v(:,q);
                % d_basis(q,:,:) = var_rec.basis.calc_basis_normal_derivatives(n1,normal,point,dbf_mag);
                d_basis1 = var_rec.basis.calc_basis_normal_derivatives(n1,normal,point,dbf_mag);
                d_basis2 = var_rec.basis.calc_basis_derivatives(n1,point,dbf_mag);
                d_basis(q,:,:) = d_basis2;
            end

            E = zeros( this.n_mtm*n_terms_local, this.n_mtm*n_terms_local );
            D = zeros( this.n_mtm*n_terms_local, this.n_mtm);
            for q = 1:n_quad
                normal = bc_fvec.v(:,q);
                % for p = 0:var_rec.degree
                % for p = var_rec.basis.degree:-1:0
                for p = 0
                    % get local projection matrix
                    % Cp = zeros(this.n_mtm,this.n_mtm);
                    % for j = 1:this.n_mtm
                    %     for i = 1:this.n_mtm
                    %         Cp(i,j) = var_rec_slip_wall_t.get_proj_rej_matrix_term(i,j,p,normal);
                    %     end
                    % end
                    Cp = -2*normal(1:this.n_mtm) * (normal(1:this.n_mtm).');
                    if (mod(p,2)~=0)
                        Cp = Cp * 0;
                    end

                    % get basis function products
                    A = zeros(n_terms_local,n_terms_local);
                    for l = 1:n_terms_local
                        for m = 1:n_terms_local
                            A(l,m) = d_basis(q,p+1,l)*d_basis(q,p+1,m);
                        end
                    end

                    % perform Kronecker product and integrate
                    E = E + bc_fquad.quad_wts(q)*kron(Cp,A);
                end

                % now for the RHS
                phi_tmp = squeeze(d_basis(q,p+1,:));
                D = D + bc_fquad.quad_wts(q)*kron(Cp,phi_tmp);
                
            end
            % E = E + kron(eye(this.n_mtm),var_rec.self_LHS);
        end
    end
    methods (Static)
        function Cijp = get_proj_rej_matrix_term(i,j,p,vec)
            ij = (i == j);
            % Cijp = ij + (-1)^p*(ij-2*vec(i)*vec(j));
            Cijp = ( (-1)^p + 1 ) * (ij + vec(i)*vec(j));
        end
    end
end