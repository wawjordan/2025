function GRID = cartesian_grid_gen(Ni,Nj,Nk)

GRID.imax = Ni;
GRID.jmax = Nj;
GRID.kmax = Nk;

[GRID.x,GRID.y,GRID.z] = ndgrid(0:Ni-1,0:Nj-1,0:Nk-1);

[GRID.xc,GRID.yc,GRID.zc] = ndgrid(1:Ni-1,1:Nj-1,1:Nk-1);
GRID.xc = GRID.xc - 0.5;
GRID.yc = GRID.yc - 0.5;
GRID.zc = GRID.zc - 0.5;

end