classdef grid_block_old < grid_block_lite_old

    properties
        nbounds   (1,1) {mustBeInteger}
        grid_vars (1,1) derived_grid_vars_t_old
        quad_ref  (1,1) reference_quad_old
        interp    (1,1) lagrange_interpolant_old
        bounds
    end

    methods
        function this = grid_block_old()
            this.interp = lagrange_interpolant_old();
        end
        function this = allocate_grid_block( this, imax, jmax, kmax )
            this = this.allocate_grid_block@grid_block_lite_old(imax,jmax,kmax);
            this.grid_vars = this.grid_vars.allocate_derived_grid( imax-1, jmax-1, kmax-1 );
            
             % will revisit this later
            this.bounds    = struct();
            for i = 1:6
                this.bounds(i).bc = constraint_container_t_old();
            end
        end
        function  this = set_bc_constraints( this, bound_num, n_var, n_dim, n_constraints, constraint_types, VAR_IDX, constraint_eval )
            this.bounds(bound_num).bc.n_constraints = n_constraints;
            this.bounds(bound_num).bc.constraint = repmat(constraint_t_old(),[n_constraints,1]);
            this.bounds(bound_num).n_var         = n_var;
            this.bounds(bound_num).n_dim         = n_dim;
            this.bounds(bound_num).constraint_types = constraint_types;
            this.bounds(bound_num).VAR_IDX = VAR_IDX;
            this.bounds(bound_num).constraint_eval = constraint_eval;
        end
        function this = fill_derived_grid_block( this, n_quad, dim )
            this = this.compute_quadrature_points( n_quad, dim );

            for k = 1:this.Ncells(3)
                for j = 1:this.Ncells(2)
                    for i = 1:this.Ncells(1)
                        this.grid_vars.volume(i,j,k) = sum( this.grid_vars.quad(i,j,k).quad_wts );
                        this.grid_vars.cell_c(1,i,j,k) = ...
                            dot( this.grid_vars.quad(i,j,k).quad_pts(1,:),...
                                 this.grid_vars.quad(i,j,k).quad_wts )    ...
                                 /this.grid_vars.volume(i,j,k);
                        this.grid_vars.cell_c(2,i,j,k) = ...
                            dot( this.grid_vars.quad(i,j,k).quad_pts(2,:),...
                                 this.grid_vars.quad(i,j,k).quad_wts )    ...
                                 /this.grid_vars.volume(i,j,k);
                        this.grid_vars.cell_c(3,i,j,k) = ...
                            dot( this.grid_vars.quad(i,j,k).quad_pts(3,:),...
                                 this.grid_vars.quad(i,j,k).quad_wts )    ...
                                 /this.grid_vars.volume(i,j,k);
                    end
                end
            end
            if dim == 2
                this.grid_vars.cell_c(3,:,:,:) = 0;
            end
            % tmp_p = { this.grid_vars.quad(:,:,:).quad_pts };
            % tmp_w = { this.grid_vars.quad(:,:,:).quad_wts };
            % this.grid_vars.volume = reshape( cellfun(@sum,tmp_w),size(this.grid_vars.volume) );
            % this.grid_vars.cell_c(1,:,:,:) = reshape( cellfun(@(A,B)A(1,:)*B,tmp_p,tmp_w), size(this.grid_vars.volume) )./this.grid_vars.volume;
            % this.grid_vars.cell_c(2,:,:,:) = reshape( cellfun(@(A,B)A(1,:)*B,tmp_p,tmp_w), size(this.grid_vars.volume) )./this.grid_vars.volume;
            % if dim == 3
            %     this.grid_vars.cell_c(3,:,:,:) = reshape( cellfun(@(A,B)A(1,:)*B,tmp_p,tmp_w), size(this.grid_vars.volume) )./this.grid_vars.volume;
            % end
            % tmp_w = { this.grid_vars.xi_face_quad(:,:,:).quad_wts };
            % this.grid_vars.xi_area = reshape( cellfun(@sum,tmp_w),size(this.grid_vars.xi_area) );
            % tmp_w = { this.grid_vars.eta_face_quad(:,:,:).quad_wts };
            % this.grid_vars.eta_area = reshape( cellfun(@sum,tmp_w),size(this.grid_vars.eta_area) );
            % if dim == 3
            %     tmp_w = { this.grid_vars.zeta_face_quad(:,:,:).quad_wts };
            %     this.grid_vars.zeta_area = reshape( cellfun(@sum,tmp_w),size(this.grid_vars.zeta_area) );
            % end
        end
        function this = compute_quadrature_points( this, n_quad, dim )
            if dim == 2
                this = compute_quadrature_points_2D( this, n_quad );
            else
                this = compute_quadrature_points_3D( this, n_quad );
            end
        end

        function this = compute_quadrature_points_2D( this, n_quad )
            this.quad_ref = reference_quad_old( [n_quad,n_quad,1] );
            quad_ref1 = quad_t_old.create_quad_ref_1D(n_quad);
            quad_ref2 = quad_t_old.create_quad_ref_2D(n_quad);
            nFine = this.Nnodes;
            nSkip = [1,1,1];
            [iCell,nCoarse] = this.get_all_coarse_indices(nFine,nSkip);
            nCellNodes = nSkip + 1;
            for m = 1:prod(nCoarse)
                iFine = this.global_to_local(iCell(m), nFine );
                iCoarse = (iFine-1)./nSkip + 1;
                [Xtmp,Ytmp,Ztmp] = this.copy_gblock_nodes(this,iFine,nSkip);
                this.grid_vars.quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = quad_t_old.map_quad_ref_2D(Xtmp,Ytmp,Ztmp,quad_ref2,this.interp);

            % xi-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(1,:,1) = true;
            this.grid_vars.xi_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_1D( ...
                            Xtmp(face_mask),...
                            Ytmp(face_mask),...
                            Ztmp(face_mask),...
                            quad_ref1,this.interp);
            this.grid_vars.xi_n( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.xi_n( iCoarse(1), iCoarse(2), iCoarse(3) ).set_face_normals_2D( Xtmp, Ytmp, Ztmp, this.quad_ref.e(5).pts, 1, this.interp );
            % xi-max faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(nCellNodes(1),:,1) = true;
            this.grid_vars.xi_face_quad( iCoarse(1)+1, iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_1D( ...
                            Xtmp(face_mask),...
                            Ytmp(face_mask),...
                            Ztmp(face_mask),...
                            quad_ref1,this.interp);
            this.grid_vars.xi_n( iCoarse(1)+1, iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.xi_n( iCoarse(1)+1, iCoarse(2), iCoarse(3) ).set_face_normals_2D( Xtmp, Ytmp, Ztmp, this.quad_ref.e(6).pts, 1, this.interp );

            % eta-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(:,1,1) = true;
            this.grid_vars.eta_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_1D( ...
                            Xtmp(face_mask),...
                            Ytmp(face_mask),...
                            Ztmp(face_mask),...
                            quad_ref1,this.interp);
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2), iCoarse(3) ).set_face_normals_2D( Xtmp, Ytmp, Ztmp, this.quad_ref.e(1).pts, 2, this.interp );
            % eta-max faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(:,nCellNodes(2),1) = true;
            this.grid_vars.eta_face_quad( iCoarse(1), iCoarse(2)+1, iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_1D( ...
                            Xtmp(face_mask),...
                            Ytmp(face_mask),...
                            Ztmp(face_mask),...
                            quad_ref1,this.interp);
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2)+1, iCoarse(3) ) = ...
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2)+1, iCoarse(3) ).set_face_normals_2D( Xtmp, Ytmp, Ztmp, this.quad_ref.e(2).pts, 2, this.interp );
            end

        end

        function this = compute_quadrature_points_3D( this, n_quad )
            this.quad_ref = reference_quad_old( repmat(n_quad,1,3) );
            quad_ref2 = quad_t_old.create_quad_ref_2D(n_quad);
            quad_ref3 = quad_t_old.create_quad_ref_3D(n_quad);
            nFine = this.Nnodes;
            nSkip = [1,1,1];
            [iCell,nCoarse] = this.get_all_coarse_indices(nFine,nSkip);
            nCellNodes = nSkip + 1;
            for m = 1:nCoarse
                iFine = this.global_to_local(iCell(m), nFine );
                iCoarse = (iFine-1)./nSkip + 1;
                [Xtmp,Ytmp,Ztmp] = this.copy_gblock_nodes(this,iFine,nSkip);
                this.grid_vars.quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = quad_t_old.map_quad_ref_3D(Xtmp,Ytmp,Ztmp,quad_ref3,this.interp);

            % xi-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(1,:,:) = true;
            this.grid_vars.xi_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_2D( ...
                            reshape(Xtmp(face_mask),nCellNodes(2),nCellNodes(3)),...
                            reshape(Ytmp(face_mask),nCellNodes(2),nCellNodes(3)),...
                            reshape(Ztmp(face_mask),nCellNodes(2),nCellNodes(3)),...
                            quad_ref2,this.interp);
            this.grid_vars.xi_n( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.xi_n( iCoarse(1), iCoarse(2), iCoarse(3) ).set_face_normals_3D( Xtmp, Ytmp, Ztmp, this.quad_ref.f(1).pts, 1, this.interp );
            % xi-max faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(nCellNodes(1),:,:) = true;
            this.grid_vars.xi_face_quad( iCoarse(1)+1, iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_2D( ...
                            reshape(Xtmp(face_mask),nCellNodes(2),nCellNodes(3)),...
                            reshape(Ytmp(face_mask),nCellNodes(2),nCellNodes(3)),...
                            reshape(Ztmp(face_mask),nCellNodes(2),nCellNodes(3)),...
                            quad_ref2,this.interp);
            this.grid_vars.xi_n( iCoarse(1)+1, iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.xi_n( iCoarse(1)+1, iCoarse(2), iCoarse(3) ).set_face_normals_3D( Xtmp, Ytmp, Ztmp, this.quad_ref.f(2).pts, 1, this.interp );

            % eta-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(:,1,:) = true;
            this.grid_vars.eta_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_2D( ...
                            reshape(Xtmp(face_mask),nCellNodes(1),nCellNodes(3)),...
                            reshape(Ytmp(face_mask),nCellNodes(1),nCellNodes(3)),...
                            reshape(Ztmp(face_mask),nCellNodes(1),nCellNodes(3)),...
                            quad_ref2,this.interp);
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2), iCoarse(3) ).set_face_normals_3D( Xtmp, Ytmp, Ztmp, this.quad_ref.f(3).pts, 2, this.interp );
            % eta-max faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(:,nCellNodes(2),:) = true;
            this.grid_vars.eta_face_quad( iCoarse(1), iCoarse(2)+1, iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_2D( ...
                            reshape(Xtmp(face_mask),nCellNodes(1),nCellNodes(3)),...
                            reshape(Ytmp(face_mask),nCellNodes(1),nCellNodes(3)),...
                            reshape(Ztmp(face_mask),nCellNodes(1),nCellNodes(3)),...
                            quad_ref2,this.interp);
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2)+1, iCoarse(3) ) = ...
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2)+1, iCoarse(3) ).set_face_normals_3D( Xtmp, Ytmp, Ztmp, this.quad_ref.f(4).pts, 2, this.interp );

            % zeta-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(:,:,1) = true;
            this.grid_vars.zeta_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t_old.map_quad_ref_2D( ...
                            reshape(Xtmp(face_mask),nCellNodes(1),nCellNodes(2)),...
                            reshape(Ytmp(face_mask),nCellNodes(1),nCellNodes(2)),...
                            reshape(Ztmp(face_mask),nCellNodes(1),nCellNodes(2)),...
                            quad_ref2,this.interp);
            this.grid_vars.zeta_n( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
            this.grid_vars.zeta_n( iCoarse(1), iCoarse(2), iCoarse(3) ).set_face_normals_3D( Xtmp, Ytmp, Ztmp, this.quad_ref.f(5).pts, 3, this.interp );
            % zeta-max faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(:,:,nCellNodes(3)) = true;
            this.grid_vars.zeta_face_quad( iCoarse(1), iCoarse(2), iCoarse(3)+1 ) = ...
                            quad_t_old.map_quad_ref_2D( ...
                            reshape(Xtmp(face_mask),nCellNodes(1),nCellNodes(2)),...
                            reshape(Ytmp(face_mask),nCellNodes(1),nCellNodes(2)),...
                            reshape(Ztmp(face_mask),nCellNodes(1),nCellNodes(2)),...
                            quad_ref2,this.interp);
            this.grid_vars.zeta_n( iCoarse(1), iCoarse(2), iCoarse(3)+1 ) = ...
            this.grid_vars.zeta_n( iCoarse(1), iCoarse(2), iCoarse(3)+1 ).set_face_normals_3D( Xtmp, Ytmp, Ztmp, this.quad_ref.f(6).pts, 3, this.interp );
            end

        end
    end % methods


    methods (Static)

        function iG = local_to_global(iSub,nSub)
            iSub = iSub(:);
            nSub = nSub(:);
            nDims = numel(iSub);
            p = 1;
            iG = 1;
            for i = 1:nDims
                iG = iG + ( iSub(i) - 1 )*p;
                p = p*nSub(i);
            end
        end
        function iSub = global_to_local(iG,nSub)
            nDims = numel(nSub);
            iSub = zeros(1,nDims);
            if (nDims==1)
                iSub(1) = iG;
                return
            end
            p = prod(nSub);
            iGtmp = iG;
            for i = nDims:-1:1
                p = fix( p/nSub(i) );
                iTmp = mod(iGtmp-1,p)+1;
                iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
                iGtmp = iTmp;
            end
        end
        function [iCoarse,nCoarse] = get_all_coarse_indices(nFine,nSkip)
            nCoarse = fix((nFine-1)./nSkip);
            ntmp = prod(nCoarse);
            iCoarse = zeros(ntmp,1);
            for i = 1:ntmp
                itmp = grid_block_old.global_to_local(i,nCoarse);
                iFine = (itmp-1).*nSkip + 1;
                iCoarse(i) = grid_block_old.local_to_global(iFine,nFine);
            end
        end
    end

end