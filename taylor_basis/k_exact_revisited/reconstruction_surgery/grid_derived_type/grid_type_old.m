classdef grid_type_old
    properties
        % nblocks is the total number of blocks in the grid
        nblocks          (1,1) {mustBeInteger}
        % total_int_cells is the total number of distinct interior cells in the grid
        total_int_cells  (1,1) {mustBeInteger}
        % total_volume is the total volume of the grid
        total_int_volume (1,1) {mustBeNumeric}
        % the block var is dimensioned to be the total number of blocks in the grid
        gblock           (:,1) grid_block_old
    end
    methods
        function this = grid_type_old(grid_struct,n_quad)
            if nargin < 2
                n_quad = 2;
            end

            if isfield(grid_struct,'nblocks')
                nblocks = grid_struct.nblocks;
            else
                nblocks = 1;
            end
            this = this.allocate_grid(nblocks);
            if isfield(grid_struct,'gblock')
                for n = 1:this.nblocks
                    imax = size(grid_struct.gblock(n).x,1);
                    jmax = size(grid_struct.gblock(n).x,2);
                    kmax = size(grid_struct.gblock(n).x,3);
                    if (kmax == 1)
                        dim = 2;
                        this.gblock(n) = this.gblock(n).allocate_grid_block(imax,jmax,2);
                        this.gblock(n).x(:,:,1) = grid_struct.gblock(n).x;
                        this.gblock(n).x(:,:,2) = grid_struct.gblock(n).x;
                        this.gblock(n).y(:,:,1) = grid_struct.gblock(n).y;
                        this.gblock(n).y(:,:,2) = grid_struct.gblock(n).y;
                        if isfield(grid_struct.gblock(n),'z')
                            this.gblock(n).z(:,:,1) = grid_struct.gblock(n).z;
                            this.gblock(n).z(:,:,2) = grid_struct.gblock(n).z + 1;
                        else
                            this.gblock(n).z(:,:,1) = 0;
                            this.gblock(n).z(:,:,2) = 1;
                        end
                    else
                        dim = 3;
                        this.gblock(n).x = grid_struct.gblock(n).x;
                        this.gblock(n).y = grid_struct.gblock(n).y;
                        this.gblock(n).z = grid_struct.gblock(n).z;
                    end
                end
            else
                imax = size(grid_struct.x,1);
                jmax = size(grid_struct.x,2);
                kmax = size(grid_struct.x,3);
                if (kmax == 1)
                    dim = 2;
                    this.gblock(1) = this.gblock(1).allocate_grid_block(imax,jmax,2);
                    this.gblock(1).x(:,:,1) = grid_struct.x;
                    this.gblock(1).x(:,:,2) = grid_struct.x;
                    this.gblock(1).y(:,:,1) = grid_struct.y;
                    this.gblock(1).y(:,:,2) = grid_struct.y;
                    if isfield(grid_struct,'z')
                        this.gblock(1).z(:,:,1) = grid_struct.z;
                        this.gblock(1).z(:,:,2) = grid_struct.z + 1;
                    else
                        this.gblock(1).z(:,:,1) = 0;
                        this.gblock(1).z(:,:,2) = 1;
                    end
                else
                    dim = 3;
                    this.gblock(1).x = grid_struct.x;
                    this.gblock(1).y = grid_struct.y;
                    this.gblock(1).z = grid_struct.z;
                end
            end
            this = this.fill_derived_grid( n_quad, dim );
        end
        function this = allocate_grid(this,nblocks)
            this.nblocks = nblocks;
            this.gblock = repmat(grid_block_old(),nblocks,1);
        end
        function this = fill_derived_grid( this, n_quad, dim )
            this.total_int_cells  = 0;
            this.total_int_volume = 0;
            for n = 1:this.nblocks
                this.gblock(n).total_cells = prod(this.gblock(n).Ncells);
                this.gblock(n) = this.gblock(n).fill_derived_grid_block(n_quad,dim);
                this.gblock(n).total_volume = sum(this.gblock(n).grid_vars.volume,"all");
                this.total_int_cells = this.total_int_cells + this.gblock(n).total_cells;
                this.total_int_volume = this.total_int_volume + this.gblock(n).total_volume;
            end
        end
    end
end