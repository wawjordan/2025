function grid = read_grd_file_to_struct( file_name )

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
    if (length(tmp) == 1)
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