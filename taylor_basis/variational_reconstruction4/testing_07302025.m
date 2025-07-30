%% Testing the var_rec_t4 type (07/30/2025)
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
agglom=false;
load_file=true;
GRID = load_gen_svf_grid_for_testing(parent_dir,agglom,load_file);

blk     = 1;
dim     = 2;
degree  = 3;
n_vars  = 4;

% [test_fun,deriv] = generate_random_poly_fun(dim,degree+1);
test_funs = cell(n_vars,1);
% for i = 1:n_vars
%     [test_funs{i},~] = generate_random_poly_fun(dim,degree+1);
% end
% for i = 1:n_vars
%     test_funs{i} = @(x,y)argEval(i,@simple_inviscid_taylor_green_vortex_2D_vel,x,y,1,0,0);
% end

gamma  = 1.4;
inputs = [2.0, 1.0, 0.8611583247416177E+05, 0.5];
% inputs = [2.0, 1.0, 0.8611583247416177E+05, 2.0;
ref_inputs = [1.0, 1.0, 347.2206293753677000 ];
test_funs{1}  = @(x,y) dens_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
test_funs{2} = @(x,y) uvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
test_funs{3} = @(x,y) vvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
test_funs{4} = @(x,y) pres_svf(x,y,0,gamma,inputs(:),ref_inputs(:));


if agglom
    idx_low  = [1,1,1];
    idx_high = [5,5,1];
else
    idx_low  = [1,1,1];
    idx_high = [10,10,1];
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

f = figure(1); set_monitor_for_figure(f,2);
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

function GRID = load_gen_svf_grid_for_testing(parent_dir,agglom,load_file)

if agglom
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_svf_agglom.mat');
else
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_svf.mat');
end

if isfile(grid_t_file)&&load_file
    load(grid_t_file,"GRID");
else
    GRID = svf_grid(11,11);
    if agglom
        GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=3,nskip=[2,2,1]);
    else
        GRID = grid_type(GRID,calc_quads=true,nquad=3);
    end
    save(grid_t_file,"GRID");
end

end

function GRID = svf_grid(N_t,N_r)

GRID.dim  = 2;
N_z = 1;
GRID.imax = N_t;
GRID.jmax = N_r;
GRID.kmax = N_z;
t_range = [pi/2, 0];
r_range = [2, 3];
z_range = [0, 0];

t = linspace(t_range(1),t_range(2),N_t);
r = linspace(r_range(1),r_range(2),N_r);
z = linspace(z_range(1),z_range(2),N_z);

[T,R,GRID.z] = ndgrid(t,r,z);

GRID.x = R.*cos(T);
GRID.y = R.*sin(T);

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