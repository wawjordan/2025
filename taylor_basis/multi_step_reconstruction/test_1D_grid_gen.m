%% testing script for "1D" grid (06/09/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = 0:10;
% x = linspace(-1,1,11);
x = rand_grd(-1,1,11,0.99);
GRID = generate_1D_grid_struct(x);
GRID = grid_type(GRID,calc_quads=true,nquad=3);




blk = 1;
dim = 1;
n1  = 1;
degree = 3;

% [test_fun,deriv] = generate_random_poly_fun(1,degree+1);
test_fun = @(x) sinc(2*x);
% test_fun = @(x) sin(2*x);
CELLS1 = set_up_cell_var_recs2(GRID,1,[1,1,1],GRID.gblock(blk).Ncells,degree,{test_fun},n1);
CELLS2 = set_up_cell_var_recs3(GRID,1,[1,1,1],GRID.gblock(blk).Ncells,degree,{test_fun},n1);
CELLS3 = CELLS2;
omega = 1.3;


stdplot(1,2);
for i = 1:100
% CELLS3 = var_rec_t2.perform_iterative_reconstruction_SOR(n1,CELLS3,omega,1);
CELLS3 = var_rec_t3.perform_iterative_reconstruction_SOR(n1,CELLS3,omega,1);
CELLS1 = var_rec_t2.perform_iterative_reconstruction_SOR(n1,CELLS1,omega,1);
% clf;
% plot_function_over_cells(test_fun,1,CELLS3(1:2),21,'k')
% plot_reconstruction_over_cells(1,CELLS2(1:2),21,'k--')
% 
% plot_reconstruction_over_cells(1,CELLS1(1:2),21,'g')
% plot_reconstruction_over_cells(1,CELLS3(1:2),21,'r--')
% legend({'exact function','exact cell averages','variational (w/o BC)', ...
%     'variational (w/ BC)'})
% drawnow
% pause(0.1)
end



stdplot(2,2);
% figure()
% exact function
plot_function_over_cells(test_fun,1,CELLS3,21,'k')

% piecewise cell avg
plot_reconstruction_over_cells(1,CELLS2,21,'k--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% old k-exact type
% stencil_size = degree+1+mod(degree,2);
% % generate stencil indices - centered in interior, biased near boundaries
% stencils = generate_simple_stencils( numel(x)-1, stencil_size );
% % generate reconstruction objects
% R = generate_reconstruction( x, degree, 1, stencils, 0 );
% % compute the reconstruction
% % retrieve cell averages
% u_avg = [CELLS1(:).coefs]; u_avg = u_avg(1,:);
% R = compute_reconstruction(R,u_avg);
% u_p1 = generate_recon_function_over_domain(R,1);
% fplot(u_p1, [x(1),x(end)],'b','LineWidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rmpath(genpath(parent_dir));
% rmpath(genpath(fullfile(parent_dir,'..')))
% addpath(genpath(fullfile(parent_dir,'..\2022_2024')))
% rmpath(genpath(fullfile(parent_dir,"'../2022_2024/SENSEI_Translation/")))
[xplt,p1,p2] = get_k_exact_rec_vals(x,degree,test_fun,21);
% rmpath(genpath(fullfile(parent_dir,'..')))
% addpath(genpath(parent_dir));

% k-exact old, multi-D
plot_reconstruction_old_over_cells(xplt,p1,'b')

% k-exact old, multi-D w/constraint
plot_reconstruction_old_over_cells(xplt,p2,'c')

% variational reconstruction without BC
plot_reconstruction_over_cells(1,CELLS1,21,'g')

% variational reconstruction with BC
plot_reconstruction_over_cells(1,CELLS3,21,'r--')


legend({'exact function','exact cell averages','k-exact (w/o BC)','k-exact (w/ BC)','variational (w/o BC)', ...
    'variational (w/ BC)'})

% [CELLS1,LHS,RHS,coefs] = var_rec_t2.perform_reconstruction(n1,CELLS1);
% figure
% plot_reconstruction_error_over_cells(test_fun,1,CELLS1,21,'k')
% plot_reconstruction_error_over_cells(test_fun,1,CELLS3,21,'r')
% plot_reconstruction_over_cells(1,CELLS2,21)

hold on;

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
[X,F] = evaluate_function_on_interp_grid(test_fun,CELLS(1).quad,npts);
plot(X{1},F{var},varargin{:})
for i = 2:numel(CELLS)
    [X,F] = evaluate_function_on_interp_grid(test_fun,CELLS(i).quad,npts);
    plot(X{1},F{var},varargin{:},'HandleVisibility','off')
end
end

function plot_reconstruction_old_over_cells(x,p,varargin)
hold on;
X = {x{1}};
F_rec= {p{1}};
plot(X{1},F_rec{1},varargin{:})
for i = 2:numel(x)
    X = {x{i}};
    F_rec= {p{i}};
    plot(X{1},F_rec{1},varargin{:},'HandleVisibility','off')
end
end

function plot_reconstruction_old_error_over_cells(x,p,test_fun,varargin)
hold on;
X = {x{1}};
F_rec= {p{1}};
plot(X{1},F_rec{1}-test_fun(X{1}),varargin{:})
for i = 2:numel(x)
    X = {x{i}};
    F_rec= {p{i}};
    plot(X{1},F_rec{1}-test_fun(X{1}),varargin{:},'HandleVisibility','off')
end
end

function plot_reconstruction_over_cells(var,CELLS,npts,varargin)
hold on;
[X,F_rec] = evaluate_reconstruction(CELLS(1),npts);
plot(X{1},F_rec{var},varargin{:})
for i = 2:numel(CELLS)
    [X,F_rec] = evaluate_reconstruction(CELLS(i),npts);
    plot(X{1},F_rec{var},varargin{:},'HandleVisibility','off')
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
function x = rand_grd(x0,xn,n,C)
% x0 - start value
% xn - stop  value
%  n - number of points
%  C - irregularity constant

% kappa = 1 + C*rand(n-1,1);
kappa = 1 + C*(1-2*rand(n-1,1));
h = (xn - x0)/sum(kappa);

x = zeros(n,1);
x(1) = x0;
x(2:n) = x0 + cumsum(kappa*h);

end