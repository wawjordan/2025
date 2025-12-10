function GRID = grid_subset(GRID,small_ind)
GRID.x = GRID.x(small_ind{:});
GRID.y = GRID.y(small_ind{:});
GRID.z = GRID.z(small_ind{:});
GRID.imax = length(small_ind{1});
GRID.jmax = length(small_ind{2});
GRID.kmax = length(small_ind{3});
end