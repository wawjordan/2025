classdef grid_vars6

    % container class for GRID type
    % has query functions to ease grabbing data
    properties
        blk      (1,1) {mustBeInteger} = 1
        n_dim    (1,1) {mustBeInteger} = 1
        lo       (:,1) {mustBeInteger}
        hi       (:,1) {mustBeInteger}
        ncells   (1,1) {mustBeInteger} = 1
        nbor_idx (:,:) {mustBeInteger}
        face_idx (:,:) {mustBeInteger}
        dir      (:,:) {mustBeInteger}
        n_int    (:,1) {mustBeInteger}
        n_bnd    (:,1) {mustBeInteger}
    end

    methods
        function this = grid_vars6( GRID, blk )
            if nargin < 1
                return
            end
            this.blk  = blk;
            if isfield( GRID, 'nblocks' ) || isprop( GRID, 'nblocks' )
                if ( blk >= 1)&&( blk <= GRID.nblocks )
                    [this.n_dim, this.lo, this.hi ] = grid_vars6.get_dim_info( GRID, blk );
                end
            else
                error("No field or property 'nblocks'")
            end

            this.ncells = prod( this.hi );
            max_nbors = 2*this.n_dim;
            this.nbor_idx = ones(this.ncells,max_nbors);
            this.face_idx = ones(this.ncells,max_nbors);
            this.dir      = ones(this.ncells,max_nbors);
            for i = 1:this.ncells
                idx = grid_vars6.global_to_local(i,this.hi);
                [this.n_int(i),this.n_bnd(i),this.dir(i,:),this.face_idx(i,:),this.nbor_idx(i,:)] = grid_vars6.get_cell_face_idx(idx, this.lo, this.hi, this.n_dim );
            end
        end

        function face_quads = get_face_quads(this,GRID,idx,interior)
            if (numel(idx)>1)
                lin_idx = grid_vars6.local_to_global_bnd(idx,this.lo,this.hi);
            else
                lin_idx = idx;
            end
            if (interior)
                face_quads = cell( this.n_int(lin_idx), 1 );
                for i = 1:this.n_int(lin_idx)
                    d = this.dir(lin_idx,i);
                    f = abs( this.face_idx(lin_idx,i) );
                    face_quads{i} = grid_vars6.get_face_quad(GRID,this.blk,d,f);
                end
            else
                face_quads = cell( this.n_bnd(lin_idx), 1 );
                for i = this.n_int(lin_idx)+1:this.n_bnd(lin_idx)
                    d = this.dir(lin_idx,i);
                    f = abs( this.face_idx(lin_idx,i) );
                    face_quads{i} = grid_vars6.get_face_quad(GRID,this.blk,d,f);
                end
            end
        end

        function face_normals = get_face_normals(this,GRID,idx,interior)
            if (numel(idx)>1)
                lin_idx = grid_vars6.local_to_global_bnd(idx,this.lo,this.hi);
            else
                lin_idx = idx;
            end
            if (interior)
                face_normals = cell( this.n_int(lin_idx), 1 );
                for i = 1:this.n_int(lin_idx)
                    d = this.dir(lin_idx,i);
                    f = this.face_idx(lin_idx,i);
                    face_normals{i} = grid_vars6.get_face_vec(GRID,this.blk,d,abs(f),sign(f));
                end
            else
                face_normals = cell( this.n_bnd(lin_idx), 1 );
                for i = this.n_int(lin_idx)+1:this.n_bnd(lin_idx)
                    d = this.dir(lin_idx,i);
                    f = this.face_idx(lin_idx,i);
                    face_normals{i} = grid_vars6.get_face_vec(GRID,this.blk,d,abs(f),sign(f));
                end
            end
        end

        function quad = get_cell_quad(this,GRID,idx)
            if (numel(idx)>1)
                lin_idx = grid_vars6.local_to_global_bnd(idx,this.lo,this.hi);
            else
                lin_idx = idx;
            end
            quad = GRID.gblock(this.blk).grid_vars.quad(lin_idx);
        end

        function h_ref = get_cell_scale(this,GRID,idx)
            if (numel(idx)==1)
                cell_idx = grid_vars6.global_to_local_bnd(idx,this.lo,this.hi);
            else
                cell_idx = idx;
            end
            nodes = grid_vars6.get_nodes(GRID,this.blk,cell_idx,this.n_dim);
            h_ref = grid_vars6.get_max_cell_extents(nodes,this.n_dim);
        end
    end

    methods (Static)
        function face_quad = get_face_quad(GRID,blk,dir,face_idx)
            quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
            face_quad = GRID.gblock(blk).grid_vars.(quad_var{dir})( face_idx );
        end
        function normals = get_face_vec(GRID,blk,dir,face_idx,factor)
            norm_var={'xi_n','eta_n','zeta_n'};
            normals = GRID.gblock(blk).grid_vars.(norm_var{dir})( face_idx );
            normals.v = normals.v * factor;
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

        function [n_int,n_bnd,dir,face_idx,nbor_idx] = get_cell_face_idx(cell_idx,bnd_min,bnd_max,dim)
            face_id = 1:2*dim;

            dir             = zeros(2*dim,1);
            linear_face_idx = zeros(2*dim,1);
            linear_nbor_idx = zeros(2*dim,1);
            for i = 1:2*dim
                [dir(i),linear_face_idx(i),linear_nbor_idx(i)] = grid_vars6.set_linear_face_idx( cell_idx, i, bnd_min(1:dim), bnd_max(1:dim) );
            end
            interior = dir>0;
            sort_idx = [face_id(interior),face_id(~interior)];
            face_idx = linear_face_idx(sort_idx);
            nbor_idx = linear_nbor_idx(sort_idx);
            dir      = abs( dir(sort_idx) );
            n_int = sum(interior);
            n_bnd = 2*dim - n_int;
        end

        function [dir,linear_face_idx,linear_nbor_idx] = set_linear_face_idx(cell_idx,face_id,bnd_min,bnd_max)

            [dir,offset,factor] = grid_vars6.face_info(face_id);
            idx_tmp = cell_idx;
            idx_tmp(dir) = idx_tmp(dir) + offset;

            nbor_idx = cell_idx;
            nbor_idx(dir) = nbor_idx(dir)+factor;

            interior = grid_vars6.in_bounds(nbor_idx,bnd_min,bnd_max);

            if (interior)
                linear_nbor_idx = grid_vars6.local_to_global_bnd(nbor_idx,bnd_min,bnd_max);
            else
                linear_nbor_idx = 0;
            end

            lo = bnd_min;
            hi = bnd_max;
            hi(dir) = hi(dir) + 1;
            linear_face_idx = factor * grid_vars6.local_to_global_bnd(idx_tmp,lo,hi);
            dir(~interior) = -dir(~interior);
        end

        function [dir,offset,factor] = face_info(face_id)
            dir    =  fix((face_id-1)/2)+1; % [1,1,2,2,3,3,...]
            offset =  mod(face_id+1,2);     % [0,1]
            factor =  2*offset-1;           % [-1, 1]
        end

        function iSub = global_to_local_bnd(iG,lo,hi)
            nSub = hi(:) - lo(:) + 1;
            iSub = grid_vars6.global_to_local(iG,nSub);
            iSub = iSub(:) + lo(:) - 1;
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

        function iG = local_to_global_bnd(iSub,lo,hi)
            iSub = iSub(:) - lo(:) + 1;
            nSub = hi(:) - lo(:) + 1;
            iG = grid_vars6.local_to_global(iSub,nSub);
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

        function val = in_bounds(idx,bnd_min,bnd_max)
            val = ( all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:)) || ...
                all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:)) );
        end

        function h_ref = get_max_cell_extents(node_list,n_dim)
            h_min = min(node_list(1:n_dim,:),[],2);
            h_max = max(node_list(1:n_dim,:),[],2);
            h_ref = 0.5 * (h_max(:) - h_min(:));
        end

        function nodes = get_nodes(GRID,blk,idx,n_dim)
            node_var={'x','y','z'};
            max_dim = 3;
            % if (n_dim > max_dim)
            %     error('GRID type not set up for more than 3 dimensions')
            % end
            n_nodes = 2^n_dim;
            nodes = zeros(max_dim,n_nodes);
            nSub = 2*ones(n_dim,1);
            for n = 1:n_nodes
                offset = grid_vars6.global_to_local(n,nSub) - 1;
                for i = 1:n_dim
                    tmp_idx = num2cell( idx(:) + offset(:) );
                    nodes(i,n) = GRID.gblock(blk).(node_var{i})( tmp_idx{:} );
                end
            end
        end
    end
end