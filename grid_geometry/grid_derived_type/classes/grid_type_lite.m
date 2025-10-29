classdef grid_type_lite
    properties
        % nblocks is the total number of blocks in the grid
        nblocks          (1,1) {mustBeInteger}
        % total_int_cells is the total number of distinct interior cells in the grid
        total_int_cells  (1,1) {mustBeInteger}
        % the block var is dimensioned to be the total number of blocks in the grid
        gblock           (:,1) grid_block_lite
    end
    methods
        function this = grid_type_lite(grid_struct)
            if nargin>0
                if isfield(grid_struct,'nblocks')
                    nb = grid_struct.nblocks;
                else
                    nb = 1;
                end
                this.nblocks = nb;
                this.gblock = repmat(grid_block_lite(),this.nblocks,1);
                this = this.copy_grid_blocks(grid_struct);
            end
        end
        function this = copy_grid_blocks(this,grid_struct)
            if isfield(grid_struct,'gblock')
                for n = 1:this.nblocks
                    imax = size(grid_struct.gblock(n).x,1);
                    jmax = size(grid_struct.gblock(n).x,2);
                    kmax = size(grid_struct.gblock(n).x,3);
                    if (kmax == 1)
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
                    this.gblock(1) = this.gblock(1).allocate_grid_block(imax,jmax,kmax);
                    this.gblock(1).x = grid_struct.x;
                    this.gblock(1).y = grid_struct.y;
                    this.gblock(1).z = grid_struct.z;
                end
            end
        end
    end
end