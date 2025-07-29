%% Testing the var_rec_t4 type (07/14/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
agglom=true;
load_file=true;
cart=true;
GRID = load_gen_grid_for_testing(parent_dir,agglom,load_file,cart);

blk     = 1;
dim     = 2;
degree  = 3;
n_vars  = 4;

% [test_fun,deriv] = generate_random_poly_fun(dim,degree+1);
test_funs = cell(n_vars,1);
% for i = 1:n_vars
%     [test_funs{i},~] = generate_random_poly_fun(dim,degree+1);
% end

for i = 1:n_vars
    test_funs{i} = @(x,y)argEval(i,@simple_inviscid_taylor_green_vortex_2D_vel,x,y,1,0,0);
end

if cart
    idx_low  = [1,1,1];
    idx_high = [10,10,1];
else
    if agglom
        % idx_low  = [30,1,1];
        % idx_high = [35,5,1];
        idx_low  = [23,1,1];
        idx_high = [42,5,1];
    else
        idx_low  = [47,1,1];
        idx_high = [84,5,1];
    end
end
SUB_GRID = GRID.subset_grid(1,idx_low,idx_high);

n1 = 1;

% CELLS_no_bc  = set_up_cell_var_recs4(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,n1,false);
% CELLS_bc     = set_up_cell_var_recs4(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,n1,true);
% CELLS2_no_bc = CELLS_no_bc;
% CELLS2_bc    = CELLS_bc;

% [CELLS_no_bc, ~,~,coefs_no_bc ] = var_rec_t4.perform_reconstruction(              n1,CELLS_no_bc );
% [CELLS2_no_bc,~,~,coefs2_no_bc] = var_rec_t4.perform_reconstruction_fully_coupled(n1,CELLS2_no_bc,false);

% [CELLS_bc, ~,~,coefs_bc ] = var_rec_t4.perform_reconstruction(              n1,CELLS_bc );
% [CELLS2_bc,~,~,coefs2_bc] = var_rec_t4.perform_reconstruction_fully_coupled(n1,CELLS2_bc,false);


CELLS  = set_up_cell_var_recs4(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,n1,false);
CELLS1 = CELLS;
CELLS2 = CELLS;
[CELLS1,~,~,coefs1] = var_rec_t4.perform_reconstruction_fully_coupled(n1,CELLS1,[]);
[CELLS2,~,~,coefs2] = var_rec_t4.perform_reconstruction_fully_coupled(n1,CELLS2,[2,3]);

%% reconstruction compare + function
% r: w/o BC
% b: w/  BC

figure(1)
var = 2;
plot_function_over_cells(test_funs{var},1,CELLS1,21,'EdgeColor','none')
% var = 3;
% plot_function_over_cells(test_funs{var},1,CELLS1,21,'EdgeColor','k')
plot_reconstruction_over_cells(var,CELLS1,21,'FaceColor','r','EdgeColor','none')
plot_reconstruction_over_cells(var,CELLS2,21,'FaceColor','b','EdgeColor','none')
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GRID = load_gen_grid_for_testing(parent_dir,agglom,load_file,cart)

grid_file = fullfile(parent_dir,'kt.grd');

if cart
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_cart.mat');
elseif agglom
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_agglom.mat');
else
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID.mat');
end

if isfile(grid_t_file)&&load_file
    load(grid_t_file,"GRID");
else
    if cart
        GRID = struct();
        GRID.twod    = true;
        GRID.nblocks = 1;
        GRID.gblock  = struct();
        GRID.gblock.imax    = 11;
        GRID.gblock.jmax    = 11;
        [GRID.gblock.x,GRID.gblock.y] = ndgrid( linspace(-1,1,GRID.gblock.imax), ...
                                                linspace(-1,1,GRID.gblock.jmax) );
        GRID = grid_type(GRID,calc_quads=true,nquad=3);
    elseif agglom
        GRID = read_grd_file_to_struct( grid_file );
        GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=3,nskip=[2,2,1]);
    else
        GRID = read_grd_file_to_struct( grid_file );
        GRID = grid_type(GRID,calc_quads=true,nquad=3);
    end
    save(grid_t_file,"GRID");
end

end


function [X,F] = evaluate_reconstruction(CELL_VAR_REC,npts,n_terms)

% get dimension
n_dim  = CELL_VAR_REC.basis.n_dim;
n_vars = CELL_VAR_REC.n_vars;

X = get_local_interp_grid(CELL_VAR_REC.quad,n_dim,npts);
F = cell(n_vars,1);
if nargin > 2
    [F{1:n_vars}] = arrayfun(@(x1,y1)CELL_VAR_REC.eval_reconstruction([x1;y1],n_terms),X{1},X{2});
else
    [F{1:n_vars}] = arrayfun(@(x1,y1)CELL_VAR_REC.eval_reconstruction([x1;y1]),X{1},X{2});
end

end

function [X,F] = evaluate_function_on_interp_grid(fun,Q,npts)
n_dim = nargin(fun);
n_var = nargout(fun);
n_var = abs(n_var);
X = get_local_interp_grid(Q,n_dim,npts);
F = cell(n_var,1);
[F{1:n_var}] = fun(X{1:n_dim});
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


function plot_function_over_cells(test_fun,var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F] = evaluate_function_on_interp_grid(test_fun,CELLS(i).quad,npts);
    surf(X{1},X{2},F{var},varargin{:})
end
end

function plot_reconstruction_over_cells(var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F_rec] = evaluate_reconstruction(CELLS(i),npts);
    surf(X{1},X{2},F_rec{var},varargin{:})
end
end
function plot_reconstruction_error_over_cells(test_fun,var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F_rec] = evaluate_reconstruction(CELLS(i),npts);
    [~,F] = evaluate_function_on_interp_grid(test_fun,CELLS(i).quad,npts);
    surf(X{1},X{2},F_rec{var}-F{var},varargin{:})
end
end

function out = argEval(n,fh,varargin)
[argsout{1:n}]=fh(varargin{:});
out = argsout{n};
end