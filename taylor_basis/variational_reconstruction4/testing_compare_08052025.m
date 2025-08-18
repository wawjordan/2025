%% Comparing var_rec_t4 with k-exact (08/05/2025)
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
% GRID = load_gen_svf_grid_for_testing(parent_dir,agglom,load_file);
GRID = load_gen_kt_grid_for_testing(parent_dir,agglom,load_file);


blk     = 1;
dim     = 2;
degree  = 5;
n_vars  = 4;

% [test_fun,deriv] = generate_random_poly_fun(dim,degree+1);
% test_funs = cell(n_vars,1);
% for i = 1:n_vars
%     [test_funs{i},~] = generate_random_poly_fun(dim,degree+1);
% end
% for i = 1:n_vars
%     test_funs{i} = @(x,y)argEval(i,@simple_inviscid_taylor_green_vortex_2D_vel,x,y,2,0,0);
% end

% gamma  = 1.4;
% inputs = [2.0, 1.0, 0.8611583247416177E+05, 0.5];
% inputs = [2.0, 1.0, 0.8611583247416177E+05, 2.0];
% ref_inputs = [1.0, 1.0, 347.2206293753677000 ];
% test_funs{1}  = @(x,y) dens_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% test_funs{2} = @(x,y) uvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% test_funs{3} = @(x,y) vvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% test_funs{4} = @(x,y) pres_svf(x,y,0,gamma,inputs(:),ref_inputs(:));


epsilon = 0.1; kappa = 0.0; tau = 0.0; alpha = 0; scale = true;
test_funs = get_airfoil_functions(epsilon,kappa,tau,alpha,scale);
test_funs = test_funs([1:3,5]);



for b = 1:6
    GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,4,2,4,[1000,1000,1000,1000],{1,2,3,5},test_funs);
end


% idx_low  = [1,1,1];
% idx_high = [10,10,1];

idx_low  = [23,1,1];
idx_high = [42,5,1];

SUB_GRID = GRID.subset_grid(1,idx_low,idx_high);

n1 = 1;

% k-exact reconstruction
CELLS1 = set_up_cell_kexact_old_recs(SUB_GRID,blk,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,false);
CELLS2 = set_up_cell_kexact_old_recs(SUB_GRID,blk,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,true);


% variational reconstruction
CELLSV1 = var_rec_t5.set_up_cell_var_recs(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,n1,false,'vector_dist');
CELLSV2 = var_rec_t5.set_up_cell_var_recs(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,n1,true,'vector_dist');

[CELLSV1,LHSV1,~,~] = var_rec_t5.perform_reconstruction_fully_coupled(n1,CELLSV1,[]);
[CELLSV2,LHSV2,~,~] = var_rec_t5.perform_reconstruction_fully_coupled(n1,CELLSV2,[]);

E_1 = get_L2_error_old(test_funs,CELLS1,CELLSV1,1:4);
E_2 = get_L2_error_old(test_funs,CELLS2,CELLSV1,1:4);

E_V1 = get_L2_error(test_funs,CELLSV1,1:4);
E_V2 = get_L2_error(test_funs,CELLSV2,1:4);

E_1_norm = squeeze(sum(E_1,[1,2]));
E_2_norm = squeeze(sum(E_2,[1,2]));
E_V1_norm = squeeze(sum(E_V1,[1,2]));
E_V2_norm = squeeze(sum(E_V2,[1,2]));

fprintf('E1:  %f  %f  %f  %f\n',E_1_norm)
fprintf('E2:  %f  %f  %f  %f\n',E_2_norm)
fprintf('EV1: %f  %f  %f  %f\n',E_V1_norm)
fprintf('EV2: %f  %f  %f  %f\n',E_V2_norm)


fprintf('condest(LHSV1) = %f\n',condest(LHSV1))
fprintf('condest(LHSV2) = %f\n',condest(LHSV2))


%% reconstruction compare + function
f1 = figure(1); clf(f1); set_monitor_for_figure(f1,2);
var = 2;

small_idx = {10:11,1:1};

plot_function_over_cells(test_funs{var},1,CELLSV1(small_idx{:}),21,'EdgeColor','none')
plot_reconstruction_over_cells_old(var,CELLSV1(small_idx{:}),CELLS1(small_idx{:}),21,'FaceColor','b','EdgeColor','none')
plot_reconstruction_over_cells_old(var,CELLSV1(small_idx{:}),CELLS2(small_idx{:}),21,'FaceColor','c','EdgeColor','none')
plot_reconstruction_over_cells(var,CELLSV1(small_idx{:}),21,'FaceColor','r','EdgeColor','none')
plot_reconstruction_over_cells(var,CELLSV2(small_idx{:}),21,'FaceColor','m','EdgeColor','none')

colorbar
axis square

hold on;

%% reconstruction error compare
f2 = figure(2); clf(f2); set_monitor_for_figure(f2,2);
var = 2;
abs_err = true;
plot_reconstruction_error_over_cells_old(abs_err,test_funs,var,CELLSV1,CELLS1,21,'FaceColor','b','EdgeColor','none','FaceAlpha',0.1)
% plot_reconstruction_error_over_cells_old(abs_err,test_funs,var,CELLSV1,CELLS2,21,'FaceColor','g','EdgeColor','none','FaceAlpha',0.1)
% plot_reconstruction_error_over_cells(abs_err,test_funs,var,CELLSV1,21,'FaceColor','r','EdgeColor','none','FaceAlpha',0.1)
plot_reconstruction_error_over_cells(abs_err,test_funs,var,CELLSV2,21,'FaceColor','c','EdgeColor','none','FaceAlpha',0.1)
colorbar
axis square
zlim([1e-6,1e-2])
set(gca,'ZScale','log')
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GRID = load_gen_kt_grid_for_testing(parent_dir,agglom,load_file)

grid_file = fullfile(parent_dir,'kt.grd');

if agglom
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_kt_agglom.mat');
else
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_kt.mat');
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


function GRID = load_gen_svf_grid_for_testing(parent_dir,agglom,load_file)

if agglom
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_svf_agglom.mat');
else
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_svf.mat');
end

if isfile(grid_t_file)&&load_file
    load(grid_t_file,"GRID");
else
    GRID = svf_grid(21,21);
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

function [E] = get_L2_error(test_fun,CELLS_VAR_REC,vars)
[E] = cell2mat(arrayfun(@(CELLS)get_cell_L2_error(test_fun,CELLS,vars,CELLS.quad),CELLS_VAR_REC,'UniformOutput',false));
end

function [E] = get_L2_error_old(test_fun,CELLS_K_REC,CELLS_VAR_REC,vars)
[E] = cell2mat(arrayfun(@(CELLS_K,CELLS_V)get_cell_L2_error_old(test_fun,CELLS_K,vars,CELLS_V.quad),CELLS_K_REC,CELLS_VAR_REC,'UniformOutput',false));
end

function [E] = get_cell_L2_error(test_fun,CELL_VAR_REC,vars,quad)
% function E = get_cell_L2_error(test_fun,CELL_VAR_REC,vars,quad)
n_var = numel(vars);
E_tmp = zeros(n_var,quad.n_quad);
for n = 1:quad.n_quad
    [E_tmp(:,n),~,~] = evaluate_reconstruction_error_at_pt(test_fun,CELL_VAR_REC,vars,quad.quad_pts(:,n));
end
E = reshape( sqrt( quad.integrate(E_tmp.^2) ),1,1,[]);
end

function [E] = get_cell_L2_error_old(test_fun,CELL_K_REC,vars,quad)
n_var = numel(vars);
E_tmp = zeros(n_var,quad.n_quad);
for n = 1:quad.n_quad
    [E_tmp(:,n),~,~] = evaluate_reconstruction_error_at_pt_old(test_fun,CELL_K_REC,vars,quad.quad_pts(:,n));
end
E = reshape( sqrt( quad.integrate(E_tmp.^2) ),1,1,[]);
end



function [E,F_rec,F_val] = evaluate_reconstruction_error_at_pt(test_funs,CELL_VAR_REC,vars,pt)

% get dimension
n_dim  = CELL_VAR_REC.basis.n_dim;
n_vars = CELL_VAR_REC.n_vars;

pt_tmp = num2cell(pt);

F = zeros(n_vars,1);
for v = 1:n_vars
    F(v) = test_funs{v}(pt_tmp{1:n_dim});
end

% index to requested variables
F_val = F(vars);
F_val = F_val(:);

F = cell(n_vars,1);
[F{1:n_vars}] = CELL_VAR_REC.eval_reconstruction(pt(1:n_dim));
F_rec = [F{vars}];
F_rec = F_rec(:);

E = F_rec - F_val;
end



function [E,F_rec,F_val] = evaluate_reconstruction_error_at_pt_old(test_funs,CELL_K_REC,vars,pt)
max_dim = 3;
% get dimension
n_dim  = CELL_K_REC.polynomial.n_dim;
n_vars = CELL_K_REC.polynomial.n_var;

pt_tmp = num2cell(pt);

F = zeros(n_vars,1);
for v = 1:n_vars
    F(v) = test_funs{v}(pt_tmp{1:n_dim});
end

% index to requested variables
F_val = F(vars);
F_val = F_val(:);

F = cell(n_vars,1);
pt_tmp = zeros(max_dim,1);
pt_tmp(1:n_dim) = pt(1:n_dim);
[F{1:n_vars}] = eval_reconstruction_old(CELL_K_REC.polynomial,pt_tmp,n_vars);
F_rec = [F{vars}];
F_rec = F_rec(:);

E = F_rec - F_val;
end


function [X,F] = evaluate_reconstruction_for_plotting(CELL_VAR_REC,npts,n_terms)

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



function [X,F] = evaluate_reconstruction_for_plotting_old(CELL_VAR_REC,CELL_K_REC,npts)

% get dimension
n_dim  = CELL_K_REC.polynomial.n_dim;
n_vars = CELL_K_REC.polynomial.n_var;

X = get_local_interp_grid(CELL_VAR_REC.quad,n_dim,npts);
F = cell(n_vars,1);
% [F{1:n_vars}] = arrayfun(@(x1,y1)CELL_K_REC.polynomial.evaluate([x1;y1;0*x1]),X{1},X{2});
[F{1:n_vars}] = arrayfun(@(x1,y1) eval_reconstruction_old(CELL_K_REC.polynomial,[x1;y1;0*x1],n_vars),X{1},X{2});
end

function varargout = eval_reconstruction_old(poly,x,n_vars)
varargout = cell(n_vars,1);
vals = poly.evaluate(x);

for v = 1:n_vars
    varargout{v} = vals(v);
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
    [X,F_rec] = evaluate_reconstruction_for_plotting(CELLS(i),npts);
    surf(X{1},X{2},F_rec{var},varargin{:})
end
end

function plot_reconstruction_error_over_cells(opt_abs,test_fun,var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F_rec] = evaluate_reconstruction_for_plotting(CELLS(i),npts);
    [~,F] = evaluate_function_on_interp_grid(test_fun{var},CELLS(i).quad,npts);
    if (opt_abs)
        surf(X{1},X{2},abs(F_rec{var}-F{1}),varargin{:})
    else
        surf(X{1},X{2},F_rec{var}-F{1},varargin{:})
    end
end
end

function plot_reconstruction_over_cells_old(var,CELLS_VAR,CELLS_K,npts,varargin)
hold on;
for i = 1:numel(CELLS_K)
    [X,F_rec] = evaluate_reconstruction_for_plotting_old( CELLS_VAR(i), CELLS_K(i), npts);
    surf(X{1},X{2},F_rec{var},varargin{:})
end
end

function plot_reconstruction_error_over_cells_old(opt_abs,test_fun,var,CELLS_VAR,CELLS_K,npts,varargin)
hold on;
for i = 1:numel(CELLS_K)
    [X,F_rec] = evaluate_reconstruction_for_plotting_old( CELLS_VAR(i), CELLS_K(i), npts );
    [~,F] = evaluate_function_on_interp_grid( test_fun{var}, CELLS_VAR(i).quad, npts );
    if (opt_abs)
        surf(X{1},X{2},abs(F_rec{var}-F{1}),varargin{:})
    else
        surf(X{1},X{2},F_rec{var}-F{1},varargin{:})
    end
end
end

function out = argEval(n,fh,varargin)
[argsout{1:n}]=fh(varargin{:});
out = argsout{n};
end