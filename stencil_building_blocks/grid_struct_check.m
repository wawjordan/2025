function GRID_OUT = grid_struct_check(grid_struct,add_centroids)
GRID_OUT = struct();
if isfield(grid_struct,'nblocks')
    nb = grid_struct.nblocks;
elseif isfield(grid_struct,'gblock')
    nb = numel(grid_struct.gblock);
else
    nb = 1;
end
GRID_OUT.nblocks = nb;
GRID_OUT = copy_grid_blocks(GRID_OUT,grid_struct,add_centroids);

GRID_OUT.twod = false;
if (any([GRID_OUT.gblock(:).twod]) )
    GRID_OUT.twod = true;
end
end

function GRID_OUT = copy_grid_blocks(GRID_OUT,grid_struct,add_centroids)
if (isfield(grid_struct,'gblock') )
    for n = 1:GRID_OUT.nblocks
        imax = size(grid_struct.gblock(n).x,1);
        jmax = size(grid_struct.gblock(n).x,2);
        kmax = size(grid_struct.gblock(n).x,3);
        if (kmax == 1)
            GRID_OUT.gblock(n) = create_gblock_struct(imax,jmax,2,true);
            GRID_OUT.gblock(n).x(:,:,1) = grid_struct.gblock(n).x;
            GRID_OUT.gblock(n).x(:,:,2) = grid_struct.gblock(n).x;
            GRID_OUT.gblock(n).y(:,:,1) = grid_struct.gblock(n).y;
            GRID_OUT.gblock(n).y(:,:,2) = grid_struct.gblock(n).y;
            if isfield(grid_struct.gblock(n),'z')
                GRID_OUT.gblock(n).z(:,:,1) = grid_struct.gblock(n).z;
                GRID_OUT.gblock(n).z(:,:,2) = grid_struct.gblock(n).z + 1;
            else
                GRID_OUT.gblock(n).z(:,:,1) = 0;
                GRID_OUT.gblock(n).z(:,:,2) = 1;
            end
            if ( add_centroids )
                [GRID_OUT.gblock(n).xc,GRID_OUT.gblock(n).yc,GRID_OUT.gblock(n).V] = quad_grid_centroids(grid_struct.gblock(n).x, grid_struct.gblock(n).y);
                GRID_OUT.gblock(n).zc       = GRID_OUT.gblock(n).z(:,:,1) + 0.5;
            end
        else
            % todo: add centroids and volume calc for 3D
            GRID_OUT.gblock(n).x = grid_struct.gblock(n).x;
            GRID_OUT.gblock(n).y = grid_struct.gblock(n).y;
            GRID_OUT.gblock(n).z = grid_struct.gblock(n).z;
        end
    end
elseif ( isfield(grid_struct,'x')&&isfield(grid_struct,'y') )
    imax = size(grid_struct.x,1);
    jmax = size(grid_struct.x,2);
    kmax = size(grid_struct.x,3);
    if (kmax == 1)
        GRID_OUT.gblock(1) = create_gblock_struct(imax,jmax,2,true);
        GRID_OUT.gblock(1).x(:,:,1) = grid_struct.x;
        GRID_OUT.gblock(1).x(:,:,2) = grid_struct.x;
        GRID_OUT.gblock(1).y(:,:,1) = grid_struct.y;
        GRID_OUT.gblock(1).y(:,:,2) = grid_struct.y;
        if isfield(grid_struct,'z')
            GRID_OUT.gblock(1).z(:,:,1) = grid_struct.z;
            GRID_OUT.gblock(1).z(:,:,2) = grid_struct.z + 1;
        else
            GRID_OUT.gblock(1).z(:,:,1) = 0;
            GRID_OUT.gblock(1).z(:,:,2) = 1;
        end
        if ( add_centroids )
            [GRID_OUT.gblock(n).xc,GRID_OUT.gblock(n).yc,GRID_OUT.gblock(n).V] = quad_grid_centroids(grid_struct.gblock(n).x, grid_struct.gblock(n).y);
            GRID_OUT.gblock(n).zc       = GRID_OUT.gblock(n).z(:,:,1) + 0.5;
        end
    else
        % todo: add centroids and volume calc for 3D
        GRID_OUT.gblock(1)   = create_gblock_struct(imax,jmax,kmax);
        GRID_OUT.gblock(1).x = grid_struct.x;
        GRID_OUT.gblock(1).y = grid_struct.y;
        GRID_OUT.gblock(1).z = grid_struct.z;
    end
else
    error('grid_struct needs to at least have "x" and "y" fields')
end
end

function gblock = create_gblock_struct(imax,jmax,kmax,twod)
gblock.imax = imax;
gblock.jmax = jmax;
gblock.kmax = kmax;
gblock.i_cells = max(imax-1,1);
gblock.j_cells = max(jmax-1,1);
gblock.k_cells = max(kmax-1,1);
gblock.Nnodes = [imax;jmax;kmax];
gblock.Ncells = max([imax-1;jmax-1;kmax-1],1);
gblock.x = zeros( imax, jmax, kmax );
gblock.y = zeros( imax, jmax, kmax );
gblock.z = zeros( imax, jmax, kmax );
if (nargin>3)
    if ( twod )
        gblock.dim = 2;
        gblock.twod = true;
    else
        gblock.dim = 3;
        gblock.twod = false;
    end
else
    if (kmax == 1 || kmax == 2)
        gblock.dim = 2;
        gblock.twod = true;
    else
        gblock.dim = 3;
        gblock.twod = false;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xc,yc,A] = quad_grid_centroids(x,y)
% order = [1,2,4,2];
% x1 = x(1:end-1,1:end-1);
% x2 = x(2:end  ,1:end-1);
% x3 = x(2:end  ,2:end  );
% x4 = x(1:end-1,2:end  );
% y1 = y(1:end-1,1:end-1);
% y2 = y(2:end  ,1:end-1);
% y3 = y(2:end  ,2:end  );
% y4 = y(1:end-1,2:end  );

[xc,yc,A] = arrayfun( @(x1,x2,x3,x4,y1,y2,y3,y4)               ...
                        quad_centroid(x1,x2,x3,x4,y1,y2,y3,y4),...
                        x(1:end-1,1:end-1), ...
                        x(2:end  ,1:end-1), ...
                        x(2:end  ,2:end  ), ...
                        x(1:end-1,2:end  ), ...
                        y(1:end-1,1:end-1), ...
                        y(2:end  ,1:end-1), ...
                        y(2:end  ,2:end  ), ...
                        y(1:end-1,2:end  ) );
end

function [xc,yc,A] = quad_centroid(x1,x2,x3,x4,y1,y2,y3,y4)
[c,A] = poly_centroid([x1,x2,x3,x4],[y1,y2,y3,y4]);
xc = c(1);
yc = c(2);
end

function [c,A] = poly_centroid(x_coords,y_coords)
x_coords = [x_coords(:);x_coords(1)];
y_coords = [y_coords(:);y_coords(1)];
N = numel(x_coords);
c = [0;0];
A = 0;
for i = 1:N-1
    det = x_coords(i)*y_coords(i+1) - x_coords(i+1)*y_coords(i);
    A    = A    + det;
    c(1) = c(1) + ( x_coords(i) + x_coords(i+1) )*det;
    c(2) = c(2) + ( y_coords(i) + y_coords(i+1) )*det;
end

c = c / (3*A);
A = A * 0.5;
end
% function A = poly_area(x_coords,y_coords)
% N = numel(x_coords);
% A = 0;
% for i = 1:N-1
%     A = A + x_coords(i)*y_coords(i+1) - x_coords(i+1)*y_coords(i);
% end
% A = 0.5 * A;
% end