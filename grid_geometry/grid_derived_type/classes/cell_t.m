classdef cell_t
    properties
        info(1,1)        block_info_t
        idx(4,1)         int32   = 0
        nbor_idx(4,6)    int32   = 0
        free(6,1)        logical = false
        bnd(6,1)         logical = false
        degree(1,1)      int32   = intmax('int32')
    end
    methods
        function this = cell_t( block, idx, block_info_list )
            if nargin ~= 3
                return
            else
                this = setup( block, idx, block_info_list );
            end
        end
        function this = setup( this, block, idx,block_info_list )
            default_offsets = this.cell_offsets();
            this.idx(1)   = block;
            this.idx(2:4) = idx;
            % get block info
            for i = 1:length(block_info_list)
                if ( block_info_list(i).block_id == this.idx(1) )
                    this.info = block_info_list(i);
                    break;
                end
            end
            % check if the current cell is on any of the block boundaries
            bnd_min = [1,1,1];
            bnd_max = this.info.Ncells;
            mask = this.face_on_a_boundary( this.idx(2:4), bnd_min, bnd_max );
            this.bnd = mask;
            for i = 1:6
                if ( mask(i) )
                    % assume out-of-bounds for now
                    this.free(i) = false;
                    % find corresponding connected boundary if it exists
                    for j = 1:length(this.info.connected_bc_list)
                        % try this BC ...
                        if ( this.info.connected_bc_list(j).my_face_id == i )
                            bnd_min1 = this.info.connected_bc_list(j).my_bnd_min;
                            bnd_max1 = this.info.connected_bc_list(j).my_bnd_max;
                            % ... if the bounds match
                            if ( this.in_bound(3,this.idx(2:4),bnd_min1,bnd_max1) )
                                % get the neighbor
                                bnd_min2 = this.info.connected_bc_list(j).nbor_bnd_min;
                                bnd_max2 = this.info.connected_bc_list(j).nbor_bnd_max;
                                this.nbor_idx(1,i) = this.info.connected_bc_list(j).nbor_block_id;
                                [ this.nbor_idx(2:4,i), ~ ] = this.get_neighbor_idx( this.idx(2:4), bnd_min1, bnd_max1, bnd_min2, bnd_max2 );
                                % check if the neighbor is in bounds
                                for k = 1:length(block_info_list)
                                    % find the corresponding neighboring
                                    % block dimensions
                                    if ( block_info_list(k).block_id == this.info.connected_bc_list(j).nbor_block_id )
                                        this.free(i) = this.in_bound( 3, this.nbor_idx(2:4,i), bnd_min, block_info_list(k).Ncells );
                                        this.bnd(i) = ~this.free(i);
                                        break;
                                    end
                                end
                            end
                        end
                    end
                else
                    % same block, simple offset
                    this.nbor_idx(:,i) = this.idx;
                    this.nbor_idx(2:4,i) = this.nbor_idx(2:4,i) + default_offsets(:,i);
                    this.free(i) = this.in_bound( 3, this.nbor_idx(2:4,i), bnd_min, bnd_max );
                end
            end
        end
    end
    methods (Static)
        function val = cell_offsets()
            val = zeros([3,6],"int32");
            direction = [-1,1];
            for d = 1:3
                for s = 1:2
                    val(d,2*(d-1)+s) = direction(s);
                end
            end
        end
        function val = in_bound(~,idx,bnd_min,bnd_max)
            val1 = all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:));
            val2 = all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:));
            val = val1||val2;
        end
        function [nbor_idx,s] = get_neighbor_idx(my_idx,idx_min1,idx_max1,idx_min2,idx_max2)
            nbor_idx = zeros(1,3);
            delta1 = (idx_max1-idx_min1);
            delta2 = (idx_max2-idx_min2);
            s1 =  1 - 2*(delta1<0);
            s2 =  1 - 2*(delta2<0);
            % can think of this as the sign for each axis across the branch cut
            s = s1.*s2;
            for i = 1:3
                offset = ( my_idx(i) - idx_min1(i) ) - s1(i);
                nbor_idx(i) = idx_min2(i) + s(i)*offset + s2(i);
            end
        end
        function mask = face_on_a_boundary(idx,bnds_min,bnds_max)
            mask = false([6,1]);
            bnds(:,1) = bnds_min;
            bnds(:,2) = bnds_max;
            for d = 1:3
                for s = 1:2
                    bnds_tmp = bnds;
                    bnds_tmp(d,:) = bnds(d,s);
                    mask(2*(d-1)+s) = cell_t.in_bound(3,idx,bnds_tmp(:,1),bnds_tmp(:,2));
                end
            end
        end
        function stencil = check_neighbors(stencil,stencil_idx,n_cells)
            % for each other cell in the stencil
            for j = 1:n_cells
                if (j==stencil_idx)
                    continue
                end
                % for each non-masked face
                for i = 1:6
                    if stencil(stencil_idx).free(i)
                        idx_tmp = stencil(stencil_idx).nbor_idx(:,i);
                        if all( stencil(j).idx == idx_tmp )
                            % update the mask on the current cell
                            stencil(stencil_idx).free(i) = false;
                            stencil(stencil_idx).degree = min( stencil(stencil_idx).degree, stencil(j).degree + 1 );
                            % and on the jth cell
                            for k = 1:6
                                if all(stencil(j).nbor_idx(:,k) == stencil(stencil_idx).idx)
                                    stencil(j).free(k) = false;
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
        function [n_faces,face_indices,face_id] = get_face_indices(dim,max_degree,stencil_cells)
            stencil_size = numel(stencil_cells);
            max_bcs = 2*dim;
            face_indices = zeros( 4, max_bcs*stencil_size );
            face_id      = zeros( 1, max_bcs*stencil_size );
            n_faces = 0;
            for s = 1:stencil_size
                if stencil_cells(s).degree > max_degree
                    continue
                end
                for n = 1:2*dim
                    if stencil_cells(s).bnd(n)
                        n_faces = n_faces + 1;
                        face_indices(1,n_faces)   = stencil_cells(s).idx(1);
                        face_indices(2:dim+1,n_faces) = stencil_cells(s).idx(2:dim+1);
                        tmp_idx = fix((n+1)/2)+1;
                        face_indices(tmp_idx,n_faces) = face_indices(tmp_idx,n_faces) + mod(n+1,2);
                        face_id(n_faces) = n;
                    end
                end
            end
        end
        function [ stencil_idx, face_idx, min_degree ] = get_grid_distance( grid, center, stencil, n_cells )
            stencil_idx = -1;
            face_idx    = -1;
            min_degree  = intmax('int32');

            degree = zeros(6*n_cells,1,'int32');
            degree(:) = intmax('int32');

            stencil_idxs = zeros(6*n_cells,1,'int32');
            min_loc_idxs = zeros(6*n_cells,1,'int32');

            distance = zeros(6*n_cells,1,'double');
            distance(:) = realmax("double");

            space_loc = zeros(1,3,'double');

            cnt = int32(0);

            % iterate through cells and calculate distances
            for j = 1:n_cells
                if ( any( stencil(j).free ) )
                    % for each neighbor
                    for i = 1:6
                        if stencil(j).free(i)

                            cnt = cnt + 1;
                            degree(cnt) = stencil(j).degree;
                            min_degree = min( degree(1:cnt) );

                            % block_id = stencil(j).idx(1);
                            % idx      = stencil(j).idx(2:4);
                            block_id = stencil(j).nbor_idx(1,i);
                            idx      = stencil(j).nbor_idx(2:4,i);
                            space_loc(:) = grid.gblock(block_id).grid_vars.cell_c( :, idx(1), idx(2), idx(3) );
                            % physical distance
                            distance(cnt) = sqrt( sum( (space_loc - center).^2 ) );
                            % fprintf('%d,%d,%f\n',j,cnt,distance(cnt));

                            stencil_idxs(cnt) = j;
                            min_loc_idxs(cnt) = i;
                        end
                    end
                end
            end

            if (cnt>0)
                if (min_degree <= 0)
                    % find the 1st minimum index
                    [~,min_idx] = min( degree(1:cnt) );

                    % Now grab any other indices that match this distance
                    mask = degree == degree(min_idx);

                    n1 = sum(mask);
                    stencil_idxs(1:n1) = stencil_idxs(mask);
                    min_loc_idxs(1:n1) = min_loc_idxs(mask);
                    distance(1:n1)   = distance(mask);
                    degree(1:n1) = degree(mask);

                    % Now find the minimum physical distance from this list
                    [~,m] = min( distance(1:n1) );
                else
                    [~,m] = min( distance(1:cnt) );
                end
                % get the offset index for that particular cell
                face_idx = min_loc_idxs(m);
                stencil_idx = stencil_idxs(m);
                min_degree = degree(m);
            end
        end
        function [ stencil, n_stencil, stencil_cells ] = build_stencil( block_id, idx, n_stencil, grid, block_info_list, balanced, print_iter )

            coord_fmt = '([%d],%3d,%3d,%3d)';
            print_fmt = sprintf('iter: %%4d, start: %s, end: %s\\n', coord_fmt,coord_fmt);

            if ( balanced )
                buffer_size    = 3*n_stencil;
                outer_iter_max = 2*buffer_size;
                inner_iter_max = 2*buffer_size;
            else
                buffer_size    = n_stencil;
                outer_iter_max = n_stencil;
                inner_iter_max = 1;
            end

            stencil_cells(buffer_size,1) = cell_t;
            center = zeros(1,3,'double');
            center(:) = grid.gblock(block_id).grid_vars.cell_c( :, idx(1), idx(2), idx(3) );

            n = 1;

            start_idx = [ block_id; idx(:) ];
            end_idx   = [ block_id; idx(:) ];

            % parent cell
            stencil_cells(n) = stencil_cells(n).setup( block_id, idx, block_info_list );
            stencil_cells(n).degree = 0;
            if ( print_iter )
                fprintf( print_fmt, n, start_idx, end_idx );
            end



            flag = false;
            outer_iter = 0;
            while ( n < n_stencil )&&( outer_iter < outer_iter_max )
                outer_iter = outer_iter + 1;
                % update masks and degree
                for i = 1:n
                    stencil_cells = cell_t.check_neighbors( stencil_cells, i, n );
                end
                % iterate through cells and calculate minimum degree
                [ stencil_idx, face_idx, min_degree ] = cell_t.get_grid_distance( grid, center, stencil_cells, n );
                cnt = 0;
                inner_iter = 0;
                current_min_degree = min_degree;
                while ( min_degree == current_min_degree )&&(inner_iter < inner_iter_max)
                    inner_iter = inner_iter + 1;
                    if ( stencil_idx == -1 )
                        warning('could not finish filling stencil, n=%d',n)
                        flag = true;
                        break;
                    end
                    if ( all(~stencil_cells(stencil_idx).free) )
                        break;
                    end

                    start_idx = stencil_cells(stencil_idx).idx;
                    end_idx   = stencil_cells(stencil_idx).nbor_idx(:,face_idx);

                    % increment the counter
                    cnt = cnt + 1;

                    % add the new cell
                    stencil_cells(n+cnt) = stencil_cells(n+cnt).setup( end_idx(1), end_idx(2:4), block_info_list );

                    % set the degree
                    stencil_cells(n+cnt).degree = min_degree + 1;

                    if ( print_iter )
                        fprintf( print_fmt, n+cnt, start_idx, end_idx );
                    end

                    % check for any edge cases
                    for i = 1:n+cnt
                        stencil_cells = cell_t.check_neighbors( stencil_cells, i, n+cnt );
                    end
                    % iterate through cells and calculate new degrees
                    [ stencil_idx, face_idx, min_degree ] = cell_t.get_grid_distance( grid, center, stencil_cells, n );
                end
                % increment the counter
                n = n + cnt;
                if flag
                    break
                end
            end
            n_stencil = n;

            stencil = zeros(4,n_stencil);
            % copy cell indices over to the stencil array
            for i = 1:n_stencil
                stencil(:,i) = stencil_cells(i).idx;
            end

        end
    end
end