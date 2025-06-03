classdef grid_type < grid_type_lite
    properties
        total_int_volume (1,1) {mustBeNumeric}
        nskip (1,3) {mustBeInteger}
    end
    methods
        function this = grid_type(grid_struct,varargin)
            this = this@grid_type_lite();
            if nargin < 1
                return
            end
            p = inputParser;
            addOptional(p,'agglomerate',false,@(x)isscalar(x)&&islogical(x))
            addOptional(p,'nskip',[1,1,1],@(x)isPositiveIntegerValuedNumeric(x)&&length(x)==3)
            addOptional(p,'mask',true,@(x)islogical(x))
            addOptional(p,'calc_quads',false,@(x)isscalar(x)&&islogical(x))
            addOptional(p,'nquad',1,@(x)isPositiveIntegerValuedNumeric(x))
            
            parse(p,varargin{:})

            if p.Results.agglomerate
                fine_grid = grid_type_lite(grid_struct);
                this = agglom_grid_blocks(this,fine_grid,p.Results.nskip);
                if p.Results.calc_quads
                    this = this.fill_derived_grid_curv( p.Results.nquad, fine_grid );
                end

                % p.Results.mask
            else
                if isfield(grid_struct,'nblocks')
                    nb = grid_struct.nblocks;
                else
                    nb = 1;
                end
                this.nblocks = nb;
                this.gblock = repmat(grid_block(),this.nblocks,1);
                this = this.copy_grid_blocks(grid_struct);
                this.nskip = [1,1,1];
                if p.Results.calc_quads
                    this = this.fill_derived_grid( p.Results.nquad );
                end
            end
        end
        function sub_grid = subset_grid(this,blk,lo,hi)
            sub_grid = grid_type();
            sub_grid.nskip   = this.nskip;
            sub_grid.nblocks = 1;
            sub_grid.gblock = this.gblock(blk).subset_grid_block(lo,hi);
            sub_grid.total_int_cells  = sub_grid.gblock(1).total_cells;
            sub_grid.total_int_volume = sub_grid.gblock(1).total_volume;
        end
        function this = agglom_grid_blocks(this,fine_grid,nskip)
            this.nblocks = fine_grid.nblocks;
            this.gblock  = repmat(grid_block(),this.nblocks,1);
            this.nskip   = nskip;
            for n = 1:this.nblocks
                imax = (fine_grid.gblock(n).imax-1)/nskip(1) + 1;
                jmax = (fine_grid.gblock(n).jmax-1)/nskip(2) + 1;
                kmax = (fine_grid.gblock(n).kmax-1)/nskip(3) + 1;
                this.gblock(n) = this.gblock(n).allocate_grid_block(imax,jmax,kmax);
                this.gblock(n).x = fine_grid.gblock(n).x(1:nskip(1):end,1:nskip(2):end,1:nskip(3):end);
                this.gblock(n).y = fine_grid.gblock(n).y(1:nskip(1):end,1:nskip(2):end,1:nskip(3):end);
                this.gblock(n).z = fine_grid.gblock(n).z(1:nskip(1):end,1:nskip(2):end,1:nskip(3):end);
            end
        end
        function this = fill_derived_grid( this, n_quad )
            this.total_int_cells  = 0;
            this.total_int_volume = 0;
            for n = 1:this.nblocks
                this.gblock(n).nskip = this.nskip;
                this.gblock(n).total_cells = prod(this.gblock(n).Ncells);
                this.gblock(n) = this.gblock(n).fill_derived_grid_block(n_quad);
                this.gblock(n).total_volume = sum(this.gblock(n).grid_vars.volume,"all");
                this.total_int_cells = this.total_int_cells + this.gblock(n).total_cells;
                this.total_int_volume = this.total_int_volume + this.gblock(n).total_volume;
            end
        end
        function this = fill_derived_grid_curv( this, n_quad, fine_grid )
            this.total_int_cells  = 0;
            this.total_int_volume = 0;
            for n = 1:this.nblocks
                this.gblock(n).nskip = this.nskip;
                this.gblock(n).total_cells = prod(this.gblock(n).Ncells);
                this.gblock(n) = this.gblock(n).fill_derived_grid_block(n_quad,fine_grid.gblock(n));
                this.gblock(n).total_volume = sum(this.gblock(n).grid_vars.volume,"all");
                this.total_int_cells = this.total_int_cells + this.gblock(n).total_cells;
                this.total_int_volume = this.total_int_volume + this.gblock(n).total_volume;
            end
        end
    end
end