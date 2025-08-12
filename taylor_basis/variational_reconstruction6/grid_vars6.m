classdef grid_vars6

    % container class for GRID type
    % has query functions to ease grabbing data
    properties
        blk    (1,1) {mustBeInteger} = 1
        n_dim  (1,1) {mustBeInteger} = 1
        lo     (:,1) {mustBeInteger}
        hi     (:,1) {mustBeInteger}
        ncells (1,1) {mustBeInteger} = 1
        idx_int (:,:) {mustBeInteger}
        n_int   (:,1)  {mustBeInteger}

    end

    methods
        function this = grid_vars6( GRID, blk )
            if nargin < 1
                return
            end
            this.blk  = blk;
            if isfield( GRID, nblocks ) || isprop( GRID, nblocks )
                if ( blk > 1)&&( blk <= GRID.nblocks )
                    [this.n_dim, this.lo, this.hi ] = grid_vars6.get_dim_info( GRID, blk );
                end
            else
                error("No field or property 'nblocks'")
            end

            this.ncells = prod( this.hi );
            max_nbors = 2*this.n_dim;
            % [idx]
            this.idx_int = ones(this.ncells,max_nbors);
            for i = 1:this.ncells
                idx = grid_vars6.global_to_local(i,this.hi);
                [nbor_idxs,cnt] = grid_vars6.get_cell_face_nbor_idx(idx, this.lo, this.hi, max_dim );
                this.n_int(i,1) = cnt;
                for n = 1:cnt
                    this.idx_int(i,n) = grid_vars6.local_to_global( nbor_idxs(:,n), this.hi );
                end
            end
        end

        function CELL = set_up_cell( this,GRID,blk,idx,degree,n_vars,dist_opt)
            [n_dim,nbor_idx] = grid_vars6.get_dim_and_nbor_info( GRID, blk );

            n_nbor = size(nbor_idx,2) - 1;
            tmp_idx = num2cell( nbor_idx(:,1) );
            volume_quad = GRID.gblock(blk).grid_vars.quad( tmp_idx{:} );
            [face_quads,face_normals,nbor_face_id,bc_face_id] = var_rec_t6_array_helper.get_face_info( GRID, ...
                blk, ...
                idx, ...
                n_dim );
            nodes = var_rec_t6_array_helper.get_nodes(GRID,blk,nbor_idx(:,1),n_dim);
            h_ref = var_rec_t6_array_helper.get_max_cell_extents(nodes,n_dim);
            basis_funs = zero_mean_basis6(n_dim,degree,h_ref,volume_quad);
            CELL = var_rec_t6_array_helper( nbor_idx, nbor_face_id(1:n_nbor), bc_face_id, ...
                basis_funs, volume_quad, face_quads, face_normals, n_vars, dist_opt );
        end

        function h_ref = get_max_cell_extents(node_list,n_dim)
            h_min = min(node_list(1:n_dim,:),[],2);
            h_max = max(node_list(1:n_dim,:),[],2);
            h_ref = 0.5 * (h_max(:) - h_min(:));
        end

        function nodes = get_nodes(GRID,blk,idx,n_dim)
            node_var={'x','y','z'};
            max_dim = 3;
            if (n_dim > max_dim)
                error('GRID type not set up for more than 3 dimensions')
            end
            n_nodes = 2^n_dim;
            nodes = zeros(max_dim,n_nodes);
            nSub = 2*ones(max_dim,1);
            for n = 1:n_nodes
                offset = zero_mean_basis6.global_to_local(n,nSub) - 1;
                for i = 1:n_dim
                    tmp_idx = num2cell( idx(:) + offset(:) );
                    nodes(i,n) = GRID.gblock(blk).(node_var{i})( tmp_idx{:} );
                end
            end
        end

        

        function [face_quads,normals,face_nbor_idx,bc_idx] = get_face_info(GRID,blk,idx,n_dim)
            quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
            norm_var={'xi_n','eta_n','zeta_n'};
            max_dim = 3;
            if (n_dim > max_dim)
                error('GRID type not set up for more than 3 dimensions')
            end
            n_faces = 2*n_dim;
            face_quads = cell(n_faces,1);
            normals    = cell(n_faces,1);
            face_nbor_idx = zeros(n_faces,1);
            bc_idx        = zeros(n_faces,1);
            cnt = 0;
            cnt_bc = 0;
            for i = 1:n_faces
                face_idx =  fix((i-1)/2)+1; % [1,1,2,2,3,3]
                offset   =  mod(i+1,2);     % [0,1]
                factor   =  2*offset-1;     % [-1, 1]
                tmp_idx  =  [idx{:}];
                tmp_idx(face_idx) = tmp_idx(face_idx)+offset;
                nbor_idx = [idx{:}];
                nbor_idx(face_idx) = nbor_idx(face_idx)+factor;
                if ( var_rec_t6_array_helper.in_bounds(           ...
                        nbor_idx(1:n_dim), ...
                        ones(1,n_dim), ...
                        GRID.gblock(blk).Ncells(1:n_dim) ) )
                    cnt = cnt + 1;
                    face_nbor_idx(cnt) = i;
                else
                    cnt_bc = cnt_bc + 1;
                    bc_idx(cnt_bc) = i;
                end
                face_quads{i} = GRID.gblock(blk).grid_vars.(quad_var{face_idx}) ...
                    ( tmp_idx(1),...
                    tmp_idx(2),...
                    tmp_idx(3) );
                normals{i}    = GRID.gblock(blk).grid_vars.(norm_var{face_idx}) ...
                    ( tmp_idx(1),...
                    tmp_idx(2),...
                    tmp_idx(3) );
                normals{i}.v = normals{i}.v * factor;
            end
            bc_idx = bc_idx(1:cnt_bc);
        end
    end

    methods (Static)
        function face_quad = get_face_quad(GRID,blk,face_id,face_idx)
            quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
            face_quad = GRID.gblock(blk).grid_vars.(quad_var{face_id})( face_idx );
        end
        function normals = get_face_normals(GRID,blk,face_id,face_idx)
            norm_var={'xi_n','eta_n','zeta_n'};
            normals = GRID.gblock(blk).grid_vars.(norm_var{face_id})( face_idx );
        end
        function [n_dim,bnd_min,bnd_max] = get_dim_info(GRID,blk)
            max_dim = 3;
            bnd_vars={'imax','jmax','kmax'};
            tmp_array = zeros(max_dim,1);
            n_dim = 0;
            for i = 1:max_dim
                if isfield( GRID.gblock(blk), bnd_vars{i} ) ...
                        || isprop( GRID.gblock(blk), bnd_vars{i} )
                    if ( GRID.gblock(blk).(bnd_vars{i}) > 2 )
                        n_dim = n_dim + 1;
                        tmp_array(n_dim) = GRID.gblock(blk).(bnd_vars{i});
                    end
                end
            end
            bnd_min = ones(n_dim,1);
            bnd_max = tmp_array(1:n_dim)-1;
        end

        function [face_nbor_idx,bc_idx] = get_face_info(GRID,blk,idx,n_dim)

            % quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
            % norm_var={'xi_n','eta_n','zeta_n'};
            n_faces = 2*n_dim;
            face_nbor_idx = zeros(n_faces,1);
            face_id       = zeros(n_faces,1);
            bc_id         = zeros(n_faces,1);
            cnt = 0;
            cnt_bc = 0;
            for i = 1:n_faces
                face_id =  fix((i-1)/2)+1; % [1,1,2,2,3,3]
                offset   =  mod(i+1,2);     % [0,1]
                factor   =  2*offset-1;     % [-1, 1]
                tmp_idx  =  [idx{:}];
                tmp_idx(face_id) = tmp_idx(face_id)+offset;
                nbor_idx = [idx{:}];
                nbor_idx(face_id) = nbor_idx(face_id)+factor;
                if ( var_rec_t6_array_helper.in_bounds(           ...
                        nbor_idx(1:n_dim), ...
                        ones(1,n_dim), ...
                        GRID.gblock(blk).Ncells(1:n_dim) ) )
                    cnt = cnt + 1;
                    face_nbor_idx(cnt) = i;
                else
                    cnt_bc = cnt_bc + 1;
                    bc_idx(cnt_bc) = i;
                end
            end
            bc_idx = bc_idx(1:cnt_bc);
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
    end
end