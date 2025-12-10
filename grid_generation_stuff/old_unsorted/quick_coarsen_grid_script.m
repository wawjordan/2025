%% Generate family of grids for Dr. Roy (07/13/2023)
clc; clear; close all;

src_file_name = 'cylinder_flare_ascii-multiblock-2D.x';
folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\TEMP_COARSEN';
jobfmt = 'cylinder_flare_ascii-multiblock-2D_%0.4dx%0.4d';

GRID = read_grd_file( fullfile(folder,src_file_name ) );
GRID = GRID.gblock;
GRID.dim = 3;
Nfine = [GRID.imax,GRID.jmax,GRID.kmax];
N_levels = 6;

for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset(GRID,{1:s:Nfine(1),1:s:Nfine(2),1:Nfine(3)});
    GRID2.dim = 3;
    filename = [folder,'\grids',sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    P3D_grid_out(GRID2,[filename,'.x'])
end

function GRID = grid_subset(GRID,small_ind)
GRID.x = GRID.x(small_ind{:});
GRID.y = GRID.y(small_ind{:});
GRID.z = GRID.z(small_ind{:});
GRID.imax = length(small_ind{1});
GRID.jmax = length(small_ind{2});
GRID.kmax = length(small_ind{3});
end

function P3D_grid_out(GRID,filename)
imax = GRID.imax;
jmax = GRID.jmax;
if GRID.kmax <= 2
    kmax = 1;
else
    kmax = GRID.kmax;
end

indices = [imax,jmax,kmax];
fid = fopen(filename,'w');
intfmt = '        %4d';
fltfmt = ' %-# 23.16E';
fprintf(fid,[intfmt,'\n'],1); % this code is just for single block grids
fprintf(fid,[repmat(intfmt,1,GRID.dim),'\n'],indices(1:GRID.dim));

count = 0;
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            fprintf(fid,fltfmt,GRID.x(i,j,k));
            count = count + 1;
            if count == 3
                fprintf(fid,'\n');
                count = 0;
            end
        end
    end
end

if GRID.dim > 1            % 2D
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                fprintf(fid,fltfmt,GRID.y(i,j,k));
                count = count + 1;
                if count == 3
                    fprintf(fid,'\n');
                    count = 0;
                end
            end
        end
    end
end

if GRID.dim > 2            % 3D
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                fprintf(fid,fltfmt,GRID.z(i,j,k));
                count = count + 1;
                if count == 3
                    fprintf(fid,'\n');
                    count = 0;
                end
            end
        end
    end
end

fclose(fid);

end


function grid = read_grd_file( file_name )

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

grid.gblock  = struct();
grid.nblocks = 1;

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.gblock.imax = str2double(tmp{1});
grid.gblock.jmax = str2double(tmp{2});
grid.gblock.kmax = 2; % 2D grid in 3D grid type

imax = grid.gblock.imax;
jmax = grid.gblock.jmax;
kmax = grid.gblock.kmax;

grid.gblock.x = zeros( imax, jmax, kmax );
grid.gblock.y = zeros( imax, jmax, kmax );
grid.gblock.z = zeros( imax, jmax, kmax );

for j = 1:jmax
    for i = 1:imax
        grid.gblock.x(i,j,1) = fscanf(fid_grid,'%f',1);
    end
end
for j = 1:jmax
    for i = 1:imax
        grid.gblock.y(i,j,1) = fscanf(fid_grid,'%f',1);
    end
end
 % copy data over to second plane
grid.gblock.x(:,:,2) = grid.gblock(1).x(:,:,1);
grid.gblock.y(:,:,2) = grid.gblock(1).y(:,:,1);
grid.gblock.z(:,:,1) = 0;
grid.gblock.z(:,:,2) = 1;

end

function grid = read_single_block_3d_grid( fid_grid )
grid = struct();
grid.nblocks = 1;

grid.gblock  = struct();

grid2d = false;

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.gblock(1).imax = str2double(tmp{1});
grid.gblock(1).jmax = str2double(tmp{2});
grid.gblock(1).kmax = str2double(tmp{3});

% if (grid.gblock(1).kmax <= 2)
%     grid2d = true;
% end

n = 1;
imax = grid.gblock(n).imax;
jmax = grid.gblock(n).jmax;
kmax = grid.gblock(n).kmax;
% if ( grid2d )
%     kmax = 1;
% end

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
% copy data over to second planes
if (grid2d)
    grid.gblock(1).kmax = 2;
    grid.gblock(1).x(:,:,2) = grid.gblock(1).x(:,:,1);
    grid.gblock(1).y(:,:,2) = grid.gblock(1).y(:,:,1); 
    grid.gblock(1).z(:,:,1) = 0;
    grid.gblock(1).z(:,:,2) = 1;
end

end

function grid = read_multi_block_2d_grid( fid_grid )
grid = struct();

line = fgetl(fid_grid);
tmp = regexp(line,'\d*','match');

grid.nblocks = str2double(tmp{1});
grid.gblock  = struct();
for n = 1:grid.nblocks
    line = fgetl(fid_grid);
    tmp = regexp(line,'\d*','match');
    
    grid.gblock(n).imax = str2double(tmp{1});
    grid.gblock(n).jmax = str2double(tmp{2});
    grid.gblock(n).kmax = 2; % 2D grid in 3D grid type
end

for n = 1:grid.nblocks
    imax = grid.gblock(n).imax;
    jmax = grid.gblock(n).jmax;
    kmax = grid.gblock(n).kmax;

    grid.gblock(n).x = zeros( imax, jmax, kmax );
    grid.gblock(n).y = zeros( imax, jmax, kmax );
    grid.gblock(n).z = zeros( imax, jmax, kmax );

    for j = 1:jmax
        for i = 1:imax
            grid.gblock(n).x(i,j,1) = fscanf(fid_grid,'%f',1);
        end
    end
    for j = 1:jmax
        for i = 1:imax
            grid.gblock(n).y(i,j,1) = fscanf(fid_grid,'%f',1);
        end
    end

    % copy data over to second planes
    grid.gblock(n).x(:,:,2) = grid.gblock(n).x(:,:,1);
    grid.gblock(n).y(:,:,2) = grid.gblock(n).y(:,:,1);
    grid.gblock(n).z(:,:,1) = 0;
    grid.gblock(n).z(:,:,2) = 1;
end

end

function grid = read_multi_block_3d_grid( fid_grid )
grid = struct();

grid2d = false;

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
%     if (grid.gblock(n).kmax <= 2)
%         grid2d = true;
%     end
end

for n = 1:grid.nblocks
    imax = grid.gblock(n).imax;
    jmax = grid.gblock(n).jmax;
    kmax = grid.gblock(n).kmax;
    

    grid.gblock(n).x = zeros( imax, jmax, kmax );
    grid.gblock(n).y = zeros( imax, jmax, kmax );
    grid.gblock(n).z = zeros( imax, jmax, kmax );

%     if (grid2d) % this breaks everything - need to figure out why
%         kmax = 1;
%     end

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
    % copy data over to second planes
    if (grid2d)
        grid.gblock(n).kmax = 2;
        grid.gblock(n).x(:,:,2) = grid.gblock(n).x(:,:,1);
        grid.gblock(n).y(:,:,2) = grid.gblock(n).y(:,:,1);
        grid.gblock(n).z(:,:,1) = 0;
        grid.gblock(n).z(:,:,2) = 1;
    end
end

end