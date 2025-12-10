function [GRIDS,F] = coarsen_refine_grid_struct( GRID_in, r_factor, refine, levels, order )
% formatting checking and regularization
GRID_in = grid_struct_check(GRID_in);

N_levels = numel(levels);

% allocate the outp grids
GRIDS = struct();
F = struct();

SZ = size_check_grid(GRID_in,r_factor,levels,refine);

tol = 1.0e-6;

if (nargin>4)
    orders = order;
else
    if (~refine)&&(all(mod(r_factor,1) < tol))
        % integer refinement
        orders = 2;
    else
        orders = 4;
    end
end


for n = 1:GRID_in.nblocks
    F.S(n) = construct_interpolant_from_grid_block( GRID_in.gblock(n), orders );
end

for j = 1:N_levels
    GRIDS(j).GRID = struct();
end

if GRID_in.twod
    for n = 1:GRID_in.nblocks
        for j = 1:N_levels
            query_vec = { linspace(-1,1,SZ(n).N(1,j)), ...
                          linspace(-1,1,SZ(n).N(2,j)) };
            GRIDS(j).GRID.gblock(n) = struct();
            GRIDS(j).GRID.gblock(n).x = F.S(n).X(query_vec);
            GRIDS(j).GRID.gblock(n).y = F.S(n).Y(query_vec);
        end
    end
else
    for n = 1:GRID_in.nblocks
        for j = 1:N_levels
            query_vec = { linspace(-1,1,SZ(n).N(1,j)), ...
                          linspace(-1,1,SZ(n).N(2,j)), ...
                          linspace(-1,1,SZ(n).N(3,j)) };
            GRIDS(j).GRID.gblock(n) = struct();
            GRIDS(j).GRID.gblock(n).x = F.S(n).X(query_vec);
            GRIDS(j).GRID.gblock(n).y = F.S(n).Y(query_vec);
            GRIDS(j).GRID.gblock(n).z = F.S(n).Z(query_vec);
        end
    end
end

for j = 1:N_levels
    GRIDS(j).GRID = grid_struct_check(GRIDS(j).GRID);
end

end

function F = construct_interpolant_from_grid_block( gblock, order )
F = struct();

F.twod = gblock.twod;
if nargin>1
    if numel(order) == 1
        orders = num2cell( order*ones(1,gblock.dim) );
    elseif nume(order) == gblock.dim
        orders = num2cell( order );
    else
        orders{1:gblock.dim} = order(1:min(gblock.dim,numel(order)));
    end
else
    orders = num2cell( 2*ones(1,gblock.dim) );
end


if F.twod
    F.Nbase = [gblock.imax,gblock.jmax];
    gridvecs = { linspace(-1,1,gblock.imax), ...
                 linspace(-1,1,gblock.jmax) };
    F.SX = spapi(orders,gridvecs,gblock.x(:,:,1));
    F.X  = @(queryvec) fnval(F.SX,queryvec);
    F.SY = spapi(orders,gridvecs,gblock.y(:,:,1));
    F.Y  = @(queryvec) fnval(F.SY,queryvec);
    F.Z  = @(queryvec) zeros(numel(queryvec{1}),numel(queryvec{2}));
else
    F.Nbase = [gblock.imax,gblock.jmax,gblock.kmax];
    gridvecs = { linspace(-1,1,gblock.imax), ...
                 linspace(-1,1,gblock.jmax), ...
                 linspace(-1,1,gblock.kmax) };
    F.SX = spapi(orders,gridvecs,gblock.x);
    F.X  = @(queryvec) fnval(F.SX,queryvec);
    F.SY = spapi(orders,gridvecs,gblock.y);
    F.Y  = @(queryvec) fnval(F.SY,queryvec);
    F.SZ = spapi(orders,gridvecs,gblock.y);
    F.Z  = @(queryvec) fnval(F.SZ,queryvec);
end
end


function SZ = size_check_grid(GRID_in,r_factor,levels,refine)
SZ = struct();

if refine
    s = +1;
else
    s = -1;
end

for n = 1:GRID_in.nblocks
    if GRID_in.gblock(n).twod
        SZ(n).N = size_check([GRID_in.gblock(n).imax,GRID_in.gblock(n).jmax],r_factor,levels,s);
    else
        SZ(n).N = size_check([GRID_in.gblock(n).imax,GRID_in.gblock(n).jmax,GRID_in.gblock(n).kmax],r_factor,levels,s);
    end
end
end

function N = size_check(Nbase,r_factor,levels,s)
sz  = numel(Nbase);
szn = numel(r_factor);
r = ones(sz,1);
if (szn==1)
    r(1:sz) = r_factor;
elseif (szn~=sz)
    r(1:min(sz,szn)) = r_factor(1:min(sz,szn));
end
N_levels = length(levels);
N = zeros(sz,N_levels);
for j = 1:N_levels
       stride = r(1:sz).^(s*levels(j));
       for i = 1:sz
           Nbase_cells = Nbase(i)-1;
           if mod( Nbase_cells*stride(i), 1 )~=0
               warning(['Warning: refinement factor %g (dimension=%d)',...
                        'at level %d results in non-integer grid size.',...
                        'Rounding to nearest.'], stride(i), i, levels(j));
               
           end
           N(i,j) = round( Nbase_cells*stride(i) ) + 1;
           if N(i,j) <= 1
               error(['Error! refinement factor %f (dimension=%d)',...
                   'at level %d is invalid'], stride(i), i, levels(j));
           end
       end
end

end

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