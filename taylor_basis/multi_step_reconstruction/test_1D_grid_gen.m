%% testing script for "1D" grid (06/09/2025)
clc; clear; close all;

% x = 0:10;
x = linspace(-1,1,11);
GRID = generate_1D_grid_struct(x);
GRID = grid_type(GRID,calc_quads=true,nquad=3);




blk = 1;
dim = 1;
n1  = 1;
degree = 3;

% [test_fun,deriv] = generate_random_poly_fun(1,degree+1);
test_fun = @(x) sinc(2*x);
CELLS2 = set_up_cell_var_recs2(GRID,1,[1,1,1],GRID.gblock(blk).Ncells,degree,{test_fun},n1);
CELLS3 = CELLS2;
omega = 1.3;

for i = 1:20
CELLS3 = var_rec_t2.perform_iterative_reconstruction_SOR(n1,CELLS3,omega,1);
clf;
plot_function_over_cells(test_fun,1,CELLS3,21,'k')
plot_reconstruction_over_cells(1,CELLS2,21,'k--')
plot_reconstruction_over_cells(1,CELLS3,21)
drawnow
pause(0.01)
end



[CELLS2,LHS,RHS,coefs] = var_rec_t2.perform_reconstruction(n1,CELLS2);

plot_reconstruction_error_over_cells(test_fun,1,CELLS2,21)
% plot_reconstruction_over_cells(1,CELLS2,21)

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
n_quad = repmat( round( Q.n_quad.^(1/ max(n_dim,2) ) ), 1, max(n_dim,2) );
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
                                                 1,...
                                                 1, quad_ref );
X = {x1,y1,z1};
end


function plot_function_over_cells(test_fun,var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F] = evaluate_function_on_interp_grid(test_fun,CELLS(i).quad,npts);
    plot(X{1},F{var},varargin{:})
end
end

function plot_reconstruction_over_cells(var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F_rec] = evaluate_reconstruction(CELLS(i),npts);
    plot(X{1},F_rec{var},varargin{:})
end
end
function plot_reconstruction_error_over_cells(test_fun,var,CELLS,npts,varargin)
hold on;
for i = 1:numel(CELLS)
    [X,F_rec] = evaluate_reconstruction(CELLS(i),npts);
    [~,F] = evaluate_function_on_interp_grid(test_fun,CELLS(i).quad,npts);
    plot(X{1},F_rec{var}-F{var},varargin{:})
end
end