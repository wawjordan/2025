function GRID = comprehensive_grid_struct(file_name)
GRID_tmp = read_grd_file_to_struct_local(file_name);
GRID = grid_struct_check_local(GRID_tmp,true);
end

function grid = read_grd_file_to_struct_local( file_name )

single_block = false;
plot2d = false;

fprintf('\n')
fprintf(' Checking grid type... ')

fid_grid = fopen(file_name,'r');

line = fgetl(fid_grid);


% Check for single block, 3D grid
tmp = regexp(line,'\d*','match');

if (length(tmp) == 3)
    plot2d       = false;
    single_block = true;
end

% Check for single block, 2D grid
if (~single_block)
    if (length(tmp) == 2)
        plot2d       = true;
        single_block = true;
    end
end

% Hopefully multiblock...
if (~single_block)
    % if (length(tmp) == 1)
    if isscalar(tmp)
        plot2d       = false;
        single_block = false;
    else
        error(['ERROR: Improperly formatted grid file! Stopping!  ',...
            'Cannot determine if grid is multiblock!  ',...
          'Check for characters or blank line(s) at beginning of file.  '])
    end
end

if (~single_block)
    plot2d = true;
    line = fgetl(fid_grid);
    tmp = regexp(line,'\d*','match');

    % Check for 3D grid
    if (length(tmp) == 3)
        plot2d = false;
    end

    % Hopefully 2D...
    if (plot2d)
        if (length(tmp) == 2)
            plot2d = true;
        else
            error(['ERROR: Improperly formatted grid file! Stopping!  ',...
                'Cannot determine if multiblock grid is 2D or 3D!'])
        end
    end
end

frewind(fid_grid);

if (plot2d)
    fprintf('%s contatins a 2D ',file_name)
    if (single_block)
        fprintf('single block grid\n');
        grid = read_single_block_2d_grid( fid_grid );
    else
        fprintf('multiblock grid\n');
        grid = read_multi_block_2d_grid( fid_grid );
    end
else
    fprintf('%s contatins a 3D ',file_name)
    if (single_block)
        fprintf('single block grid\n');
        grid = read_single_block_3d_grid( fid_grid );
    else
        fprintf('multiblock grid\n');
        grid = read_multi_block_3d_grid( fid_grid );
    end
end
fclose(fid_grid);

end

function grid = read_single_block_2d_grid( fid_grid )
grid = struct();

grid.twod = true;

grid.gblock  = struct();
grid.nblocks = 1;

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.gblock.imax = str2double(tmp{1});
grid.gblock.jmax = str2double(tmp{2});

imax = grid.gblock.imax;
jmax = grid.gblock.jmax;

grid.gblock.x = zeros( imax, jmax );
grid.gblock.y = zeros( imax, jmax );

for j = 1:jmax
    for i = 1:imax
        grid.gblock.x(i,j) = fscanf(fid_grid,'%f',1);
    end
end
for j = 1:jmax
    for i = 1:imax
        grid.gblock.y(i,j) = fscanf(fid_grid,'%f',1);
    end
end

end

function grid = read_single_block_3d_grid( fid_grid )
grid = struct();

grid.twod = false;

grid.nblocks = 1;

grid.gblock  = struct();

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.gblock(1).imax = str2double(tmp{1});
grid.gblock(1).jmax = str2double(tmp{2});
grid.gblock(1).kmax = str2double(tmp{3});

n = 1;
imax = grid.gblock(n).imax;
jmax = grid.gblock(n).jmax;
kmax = grid.gblock(n).kmax;

grid.gblock(n).x = zeros( imax, jmax, kmax );
grid.gblock(n).y = zeros( imax, jmax, kmax );
grid.gblock(n).z = zeros( imax, jmax, kmax );

for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            grid.gblock(1).x(i,j,k) = fscanf(fid_grid,'%f',1);
        end
    end
end
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            grid.gblock(1).y(i,j,k) = fscanf(fid_grid,'%f',1);
        end
    end
end
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            grid.gblock(1).z(i,j,k) = fscanf(fid_grid,'%f',1);
        end
    end
end

end

function grid = read_multi_block_2d_grid( fid_grid )
grid = struct();

grid.twod = true;

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.nblocks = str2double(tmp{1});
grid.gblock  = struct();
for n = 1:grid.nblocks
    line = fgetl(fid_grid);
    tmp = regexp(line,'\d*','match');
    
    grid.gblock(n).imax = str2double(tmp{1});
    grid.gblock(n).jmax = str2double(tmp{2});
end

for n = 1:grid.nblocks
    imax = grid.gblock(n).imax;
    jmax = grid.gblock(n).jmax;

    grid.gblock(n).x = zeros( imax, jmax );
    grid.gblock(n).y = zeros( imax, jmax );

    for j = 1:jmax
        for i = 1:imax
            grid.gblock(n).x(i,j) = fscanf(fid_grid,'%f',1);
        end
    end
    for j = 1:jmax
        for i = 1:imax
            grid.gblock(n).y(i,j) = fscanf(fid_grid,'%f',1);
        end
    end
end

end

function grid = read_multi_block_3d_grid( fid_grid )
grid = struct();

grid.twod = false;

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.nblocks = str2double(tmp{1});
grid.gblock  = struct();

for n = 1:grid.nblocks
    line = fgetl(fid_grid);
    tmp = regexp(line,'\d*','match');
    
    grid.gblock(n).imax = str2double(tmp{1});
    grid.gblock(n).jmax = str2double(tmp{2});
    grid.gblock(n).kmax = str2double(tmp{3});
end

for n = 1:grid.nblocks
    imax = grid.gblock(n).imax;
    jmax = grid.gblock(n).jmax;
    kmax = grid.gblock(n).kmax;
    

    grid.gblock(n).x = zeros( imax, jmax, kmax );
    grid.gblock(n).y = zeros( imax, jmax, kmax );
    grid.gblock(n).z = zeros( imax, jmax, kmax );

    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                grid.gblock(n).x(i,j,k) = fscanf(fid_grid,'%f',1);
            end
        end
    end
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                grid.gblock(n).y(i,j,k) = fscanf(fid_grid,'%f',1);
            end
        end
    end
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                grid.gblock(n).z(i,j,k) = fscanf(fid_grid,'%f',1);
            end
        end
    end
end

end

function GRID_OUT = grid_struct_check_local(grid_struct,add_centroids)
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
            [GRID_OUT.gblock(1).xc,GRID_OUT.gblock(1).yc,GRID_OUT.gblock(1).V] = quad_grid_centroids(grid_struct.gblock(1).x, grid_struct.gblock(1).y);
            GRID_OUT.gblock(1).zc = GRID_OUT.gblock(1).z(:,:,1) + 0.5;
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
[c,A,~] = poly_centroid([x1,x2,x3,x4],[y1,y2,y3,y4]);
xc = c(1);
yc = c(2);
end

function [c,A,winding] = poly_centroid(x_coords,y_coords)
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

winding = A>0;
A = abs(A) * 0.5;
end
% function A = poly_area(x_coords,y_coords)
% N = numel(x_coords);
% A = 0;
% for i = 1:N-1
%     A = A + x_coords(i)*y_coords(i+1) - x_coords(i+1)*y_coords(i);
% end
% A = 0.5 * A;
% end