function GRID_OUT = grid_struct_check(grid_struct)
GRID_OUT = struct();
if isfield(grid_struct,'nblocks')
    nb = grid_struct.nblocks;
elseif isfield(grid_struct,'gblock')
    nb = numel(grid_struct.gblock);
else
    nb = 1;
end
GRID_OUT.nblocks = nb;
GRID_OUT = copy_grid_blocks(GRID_OUT,grid_struct);

GRID_OUT.twod = false;
if (any([GRID_OUT.gblock(:).twod]) )
    GRID_OUT.twod = true;
end
end

function GRID_OUT = copy_grid_blocks(GRID_OUT,grid_struct)
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
        else
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
    else
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