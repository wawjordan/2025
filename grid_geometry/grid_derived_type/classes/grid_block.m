classdef grid_block < grid_block_lite

    properties
        nbounds   (1,1) {mustBeInteger}
        grid_vars (1,1) derived_grid_vars_t
        quad_ref  (1,1) reference_quad
        interp    (1,1) lagrange_interpolant
        nskip     (1,3) {mustBeInteger} = [1,1,1]
        dim       (1,1) {mustBeInteger}
        bounds
    end

    methods
        function this = grid_block(nskip)
            this.interp = lagrange_interpolant();
            if nargin>0
                this.nskip=nskip;
            end
        end
        function sub_grid = subset_grid_block(this,lo,hi)
            sub_grid = grid_block(this.nskip);
            sub_grid.imax = hi(1) - lo(1) + 2;
            sub_grid.jmax = hi(2) - lo(2) + 2;
            sub_grid.kmax = hi(3) - lo(3) + 2;
            sub_grid.i_cells = sub_grid.imax-1;
            sub_grid.j_cells = sub_grid.jmax-1;
            sub_grid.k_cells = sub_grid.kmax-1;
            sub_grid.Nnodes = [sub_grid.imax;sub_grid.jmax;sub_grid.kmax];
            sub_grid.Ncells = [sub_grid.imax-1;sub_grid.jmax-1;sub_grid.kmax-1];

            sub_grid.total_cells = prod(sub_grid.Ncells);

            sub_grid.x  = this.x(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1);
            sub_grid.y  = this.y(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1);
            sub_grid.z  = this.z(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1);
            sub_grid.allocated = true;
            sub_grid.dim = this.dim;
            sub_grid.nbounds  = this.nbounds; %%%%%
            sub_grid.bounds   = this.bounds;
            sub_grid.quad_ref = this.quad_ref;
            sub_grid.block_id = this.block_id;

            sub_grid.grid_vars = sub_grid.grid_vars.allocate_derived_grid( sub_grid.imax-1,sub_grid.jmax-1,sub_grid.kmax-1 );
            sub_grid.grid_vars.cell_c    = this.grid_vars.cell_c( :,lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)  );
            sub_grid.grid_vars.volume    = this.grid_vars.volume(   lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)  );
            sub_grid.grid_vars.quad      = this.grid_vars.quad(     lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)  );
            sub_grid.grid_vars.xi_area   = this.grid_vars.xi_area(  lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  );
            sub_grid.grid_vars.eta_area  = this.grid_vars.eta_area( lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  );
            sub_grid.grid_vars.zeta_area = this.grid_vars.zeta_area(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1);
            sub_grid.grid_vars.xi_n   = this.grid_vars.xi_n(  lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  );
            sub_grid.grid_vars.eta_n  = this.grid_vars.eta_n( lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  );
            sub_grid.grid_vars.zeta_n = this.grid_vars.zeta_n(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1);
            sub_grid.grid_vars.xi_face_quad   = this.grid_vars.xi_face_quad(  lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  );
            sub_grid.grid_vars.eta_face_quad  = this.grid_vars.eta_face_quad( lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  );
            sub_grid.grid_vars.zeta_face_quad = this.grid_vars.zeta_face_quad(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1);

            sub_grid.total_volume = sum(sub_grid.grid_vars.volume,"all");
            
        end
        function this = allocate_grid_block( this, imax, jmax, kmax )
            this = this.allocate_grid_block@grid_block_lite(imax,jmax,kmax);
            if (kmax == 1 || kmax == 2)
                this.dim = 2;
            else
                this.dim = 3;
            end
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
        function this = fill_derived_grid_block( this, n_quad, varargin )
            p = inputParser;
            p.addOptional('fine_grid_block',grid_block_lite())
            parse(p,varargin{:})
            if p.Results.fine_grid_block.allocated
                this = this.compute_quadrature_points( n_quad, p.Results.fine_grid_block );
            else
                this = this.compute_quadrature_points( n_quad );
            end

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
            if this.dim == 2
                this.grid_vars.cell_c(3,:,:,:) = 0;
            end
        end
        function this = compute_quadrature_points( this, n_quad, fine_grid_block )
            if nargin>2
                if this.dim == 2
                    this = compute_quadrature_points_2D( this, n_quad, fine_grid_block );
                else
                    this = compute_quadrature_points_3D( this, n_quad, fine_grid_block );
                end
            else
                if this.dim == 2
                    this = compute_quadrature_points_2D( this, n_quad );
                else
                    this = compute_quadrature_points_3D( this, n_quad );
                end
            end
        end

        function this = compute_quadrature_points_2D( this, n_quad, fine_grid_block )
            this.quad_ref = reference_quad( [n_quad,n_quad,1] );
            quad_ref1 = quad_t.create_quad_ref_1D(n_quad);
            quad_ref2 = quad_t.create_quad_ref_2D(n_quad);

            block = this;
            if nargin>2
                block = fine_grid_block;
            end
            nFine = block.Nnodes;
            
            nSkip = this.nskip;
            [iCell,nCoarse] = this.get_all_coarse_indices(nFine,nSkip);
            nCellNodes = nSkip + 1;

            for m = 1:prod(nCoarse)
                iFine = this.global_to_local(iCell(m), nFine );
                iCoarse = (iFine-1)./nSkip + 1;
                [Xtmp,Ytmp,Ztmp] = this.copy_gblock_nodes(block,iFine,nSkip);
                this.grid_vars.quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = quad_t.map_quad_ref_2D(Xtmp,Ytmp,Ztmp,quad_ref2,this.interp);

            % xi-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(1,:,1) = true;
            this.grid_vars.xi_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t.map_quad_ref_1D( ...
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
                            quad_t.map_quad_ref_1D( ...
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
                            quad_t.map_quad_ref_1D( ...
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
                            quad_t.map_quad_ref_1D( ...
                            Xtmp(face_mask),...
                            Ytmp(face_mask),...
                            Ztmp(face_mask),...
                            quad_ref1,this.interp);
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2)+1, iCoarse(3) ) = ...
            this.grid_vars.eta_n( iCoarse(1), iCoarse(2)+1, iCoarse(3) ).set_face_normals_2D( Xtmp, Ytmp, Ztmp, this.quad_ref.e(2).pts, 2, this.interp );
            end

        end

        function this = compute_quadrature_points_3D( this, n_quad, fine_grid_block )
            this.quad_ref = reference_quad( repmat(n_quad,1,3) );
            quad_ref2 = quad_t.create_quad_ref_2D(n_quad);
            quad_ref3 = quad_t.create_quad_ref_3D(n_quad);

            block = this;
            if nargin>2
                block = fine_grid_block;
            end
            nFine = block.Nnodes;
            
            nSkip = this.nskip;
            [iCell,nCoarse] = this.get_all_coarse_indices(nFine,nSkip);
            nCellNodes = nSkip + 1;
            
            for m = 1:prod(nCoarse)
                iFine = this.global_to_local(iCell(m), nFine );
                iCoarse = (iFine-1)./nSkip + 1;
                [Xtmp,Ytmp,Ztmp] = this.copy_gblock_nodes(block,iFine,nSkip);
                this.grid_vars.quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = quad_t.map_quad_ref_3D(Xtmp,Ytmp,Ztmp,quad_ref3,this.interp);

            % xi-min faces
            face_mask = false(nCellNodes(1),nCellNodes(2),nCellNodes(3));
            face_mask(1,:,:) = true;
            this.grid_vars.xi_face_quad( iCoarse(1), iCoarse(2), iCoarse(3) ) = ...
                            quad_t.map_quad_ref_2D( ...
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
                            quad_t.map_quad_ref_2D( ...
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
                            quad_t.map_quad_ref_2D( ...
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
                            quad_t.map_quad_ref_2D( ...
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
                            quad_t.map_quad_ref_2D( ...
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
                            quad_t.map_quad_ref_2D( ...
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
                itmp = grid_block.global_to_local(i,nCoarse);
                iFine = (itmp-1).*nSkip + 1;
                iCoarse(i) = grid_block.local_to_global(iFine,nFine);
            end
        end
    end

end