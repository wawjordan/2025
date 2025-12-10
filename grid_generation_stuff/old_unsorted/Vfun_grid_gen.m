function GRID = Vfun_grid_gen(Vfun,Ni,Nj,Nk,ranges)

GRID.imax = Ni;
GRID.jmax = Nj;
GRID.kmax = Nk;

[X1,X2,X3] = ndgrid(linspace(ranges(1,1),ranges(1,2),Ni),...
                    linspace(ranges(2,1),ranges(2,2),Nj),...
                    linspace(ranges(3,1),ranges(3,2),Nk) );

GRID.x = Vfun.X(X1,X2,X3);
GRID.y = Vfun.Y(X1,X2,X3);
GRID.z = Vfun.Z(X1,X2,X3);

end