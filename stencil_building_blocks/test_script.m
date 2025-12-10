%% script for generating a basic stencil (2D)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = 'stencil_building_blocks';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;


%% load the grid
% when run for the first time, it will take a while to calculate everything in
% the grid derived type. Therefore, if you are just exploring the same grid,
% you can set load_file=true to load the corresponding .mat file instead of
% re-calculating everything
load_file = true;

% whether or not to agglomerate the mesh. Currently hardcoded with
% nskip=[2,2,1]
agglom   = true;

grid_file = 'kt.grd';

% generate (or load) a grid class for the specified geometry
GRID = load_gen_grid_for_testing(parent_dir,grid_file,agglom,load_file);

%% parameters for stencil growth
% block selection
blk = 1;

% index selection (center idx of your stencil)
idx = [16,1,1];

% number of requested cells in the stencil
n_stencil = 8;

% whether or not the stencil should be 'balanced'
% i.e. all cells of a given topological distance away are selected
balanced = true;

% whether or not to print the starting and ending idxs during stencil growth
print_iter = true;

% block information (needed to inform stencil growth near boundaries
block_info_list(1)           = block_info_t();
block_info_list(1).block_id  = 1;
block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;

% this is specifically for the airfoil (kt.grd) with agglom=true
% comment it out if using a different grid
block_info_list(1).connected_bc_list(1)                 = connected_bc_info_t();
block_info_list(1).connected_bc_list(1).my_block_id     = 1;
block_info_list(1).connected_bc_list(1).my_face_id      = 3;
block_info_list(1).connected_bc_list(1).nbor_face_id    = 3;
block_info_list(1).connected_bc_list(1).my_bnd_min(:)   = [  1, 1, 1 ];
block_info_list(1).connected_bc_list(1).my_bnd_max(:)   = [ 24, 1, 1 ];
block_info_list(1).connected_bc_list(1).nbor_block_id   = 1;
block_info_list(1).connected_bc_list(1).nbor_bnd_min(:) = [64, 1, 1 ];
block_info_list(1).connected_bc_list(1).nbor_bnd_max(:) = [ 41, 1, 1 ];

block_info_list(1).connected_bc_list(2)                 = connected_bc_info_t();
block_info_list(1).connected_bc_list(2).my_block_id     = 1;
block_info_list(1).connected_bc_list(2).my_face_id      = 3;
block_info_list(1).connected_bc_list(2).nbor_face_id    = 3;
block_info_list(1).connected_bc_list(2).my_bnd_min(:)   = [ 41, 1, 1 ];
block_info_list(1).connected_bc_list(2).my_bnd_max(:)   = [64, 1, 1 ];
block_info_list(1).connected_bc_list(2).nbor_block_id   = 1;
block_info_list(1).connected_bc_list(2).nbor_bnd_min(:) = [ 24, 1, 1 ];
block_info_list(1).connected_bc_list(2).nbor_bnd_max(:) = [  1, 1, 1 ];

%% build the stencil
[ stencil_idx, n_stencil_out, stencil_cells ] = cell_t.build_stencil( blk, idx, n_stencil, GRID, block_info_list, balanced, print_iter );

%% plot and label the stencil
plot_stencil_cells_2D(GRID,stencil_idx,'k')
label_stencil_cells(GRID,stencil_idx,true)
label_stencil_cells(GRID,stencil_idx,false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GRID = load_gen_grid_for_testing(parent_dir,grid_file,agglom,load_file)

grid_t_file = [extractBefore(grid_file,'.'),'_GRID'];

grid_file = fullfile(parent_dir,grid_file);

if agglom
    grid_t_file = fullfile(parent_dir,[grid_t_file,'_agglom.mat']);
else
    grid_t_file = fullfile(parent_dir,[grid_t_file,'.mat']);
end

if isfile(grid_t_file)&&load_file
    load(grid_t_file,"GRID");
else
    GRID = read_grd_file_to_struct(grid_file);
    if agglom
        GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=3,nskip=[2,2,1]);
    else
        GRID = grid_type(GRID,calc_quads=true,nquad=3);
    end
    save(grid_t_file,"GRID");
end

end

function grid_in_stencil = get_stencil_grid(Npts,stencil,grid)
n_stencil = size(stencil,2);
grid_in_stencil = cell(n_stencil,3);
for si = 1:n_stencil
    b   = stencil(1,si);
    idx = stencil(2:4,si);
    idx = idx(:).';
    X = get_local_interp_grid(grid.gblock(b).grid_vars.quad(idx(1),idx(2),idx(3)), grid.gblock(b).dim,Npts);
    grid_in_stencil{si,1} = X{1};
    grid_in_stencil{si,2} = X{2};
    grid_in_stencil{si,3} = X{3};
end
end

function X = get_local_interp_grid(Q,n_dim,npts)
% get number of quadrature points (per dimension)
n_quad = repmat( round( Q.n_quad.^(1/n_dim) ), 1, n_dim );
quad_ref = quad_t.create_quad_ref_2D( n_quad(1) );


n_plot_pts = ones(3,1);
n_plot_pts(1:n_dim) = npts;

xtmp = reshape(Q.quad_pts(1,:),n_quad);
ytmp = reshape(Q.quad_pts(2,:),n_quad);
ztmp = reshape(Q.quad_pts(2,:),n_quad);

L = lagrange_interpolant(n_quad(1) + 1);
% extrapolate the parametric space from the innermost quadrature points
% Note: needs to be a consistent order of quadrature for the geometry representation
[x1,y1,z1] = L.map_grid_from_quad(xtmp,ytmp,ztmp,n_plot_pts(1),...
                                                 n_plot_pts(2),...
                                                 n_plot_pts(3), quad_ref );
X = {x1,y1,z1};
end

function plot_stencil_cells_2D(GRID,stencil,varargin)

k_plane = 1;
Npts = 2;
if any(GRID.nskip>1)
    Npts = 21;
end
grid_in_stencil = get_stencil_grid(Npts,stencil,GRID);
n_stencil = size(stencil,2);

hold on;
for n = 1:n_stencil
    X = grid_in_stencil{n,1}(:,:,k_plane);
    Y = grid_in_stencil{n,2}(:,:,k_plane);
    Z = grid_in_stencil{n,3}(:,:,k_plane);
    plot_2D_mesh_edge(X,Y,Z,varargin{:});
end

end

function label_stencil_cells(GRID,stencil,print_idx)
n_stencil = size(stencil,2);
if (print_idx)
    for n = 1:n_stencil
        b = stencil(1,n);
        idx = stencil(2:4,n);
        label_cell(GRID.gblock(b).grid_vars.cell_c(:,idx(1),idx(2),idx(3)),'%d:[%d,%d,%d]',[n;idx(:)],'Color','k')
    end
else
    for n = 1:n_stencil
        b = stencil(1,n);
        idx = stencil(2:4,n);
        label_cell(GRID.gblock(b).grid_vars.cell_c(:,idx(1),idx(2),idx(3)),'%d',n,'Color','r')
    end
end
    
end



function p = plot_2D_mesh_edge(X,Y,Z,varargin)
hold on
X1 = X(1:end,1);
X2 = X(end,1:end);
X3 = X(end:-1:1,end);
X4 = X(1,end:-1:1);
Xp = [X1(:);X2(:);X3(:);X4(:)];

Y1 = Y(1:end,1);
Y2 = Y(end,1:end);
Y3 = Y(end:-1:1,end);
Y4 = Y(1,end:-1:1);
Yp = [Y1(:);Y2(:);Y3(:);Y4(:)];

Z1 = Z(1:end,1);
Z2 = Z(end,1:end);
Z3 = Z(end:-1:1,end);
Z4 = Z(1,end:-1:1);
Zp = [Z1(:);Z2(:);Z3(:);Z4(:)];
p = plot3( Xp, Yp, Zp, varargin{:} );
end

function label_cell(X,fmt,val,varargin)
text(X(1),X(2),X(3),sprintf(fmt,val),varargin{:})
end