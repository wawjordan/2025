function GRID = cylinder_grid_gen(N_r,N_t,N_z,r_range,t_range,z_range,spacing)
% spacing is a cell array of spacing functions for the grid

GRID.imax = N_r;
GRID.jmax = N_t;
GRID.kmax = N_z;

r = spacing{1}(r_range(1),r_range(2),N_r);
t = spacing{2}(t_range(1),t_range(2),N_t);
z = spacing{3}(z_range(1),z_range(2),N_z);

[R,T,GRID.z] = ndgrid(r,t,z);

GRID.x = R.*cos(T);
GRID.y = R.*sin(T);

end