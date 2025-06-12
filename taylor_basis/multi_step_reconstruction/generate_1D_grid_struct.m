function GRID = generate_1D_grid_struct(x)

GRID = struct();
GRID.twod    = true;
GRID.nblocks = 1;
GRID.gblock  = struct();
GRID.gblock.imax = numel(x);
GRID.gblock.jmax = 2;
[GRID.gblock.x,GRID.gblock.y] = ndgrid( x, [0,1] );

end