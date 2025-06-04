%% Testing the taylor_shape_functions type (04/30/2025)
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
cart=false;
GRID = load_gen_grid_for_testing(parent_dir,agglom,load_file,cart);

blk     = 1;
dim     = 2;
degree  = 4;
n_vars  = 1;
[test_fun,test_fun_grad,test_fun_hess] = generate_test_function();

if cart
    idx_low  = [1,1,1];
    idx_high = [5,5,1];
else
    if agglom
        idx_low  = [23,1,1];
        idx_high = [42,5,1];
    else
        idx_low  = [47,1,1];
        idx_high = [84,5,1];
    end
end
SUB_GRID = GRID.subset_grid(1,idx_low,idx_high);

CELLS = set_up_cells(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,n_vars);
CELLS = arrayfun(@(CELLS)CELLS.set_cell_avg({test_fun}),CELLS);

n1 = 1;
CELLS2 = set_up_cell_var_recs2(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,{test_fun},n1);
[CELLS2,LHS,RHS,coefs] = var_rec_t2.perform_reconstruction(n1,CELLS2);
% LHS = full(LHS);
% LHS(abs(LHS)<1e-6)=0;

plot_reconstruction_error_over_cells(test_fun,1,CELLS2,21,'EdgeColor','none')
view(12,14)

T1 = CELLS(1).taylor;
Q1 = CELLS(1).quad;

T2 = CELLS2(1).basis;
Q2 = CELLS2(1).quad;

coefs1 = CELLS2(1).coefs;
coefs2 = get_exact_quadratic_reconstruction_in_cell(T1,Q1,test_fun,test_fun_grad,test_fun_hess);

[X1,F1] = evaluate_function_on_interp_grid(test_fun,CELLS2(1).quad,21);
[~,F_rec1] = evaluate_reconstruction(CELLS2(1),21);

[X2,F2] = evaluate_function_on_interp_grid(test_fun,CELLS2(2).quad,21);
[~,F_rec2] = evaluate_reconstruction(CELLS2(2),21);

% CELLS2(1).coefs = coefs2(1:CELLS2(1).basis.n_terms);
% [X2,F_rec2] = evaluate_reconstruction(CELLS2(1),21);

hold on;
surf(X1{1},X1{2},F1{1},'EdgeColor','k')
view(-14,12)
surf(X2{1},X2{2},F2{1},'EdgeColor','k')
surf(X1{1},X1{2},F_rec1{1},'FaceColor','b','EdgeColor','k')
surf(X2{1},X2{2},F_rec2{1},'FaceColor','r','EdgeColor','k')    
point = CELLS(1).fquad(1).quad_pts(:,1);

% l = 5;
% [tmp1,delta1] = T1.eval_derivative( 5, point, T1.exponents(:,2) );
% [tmp2,delta2] = T2.deval( 5, point, T2.terms(2).exponents );



figure()
C = colororder(gcf);
for n = 1:T1.n_terms
    [X1,F1,~] = evaluate_taylor_basis_1(T1,n,Q1,21);
    [X2,F2,~] = evaluate_taylor_basis_2(T2,n,Q2,21);
    hold on;
    surf(X1{1},X1{2},F1,'FaceColor','b','EdgeColor','k')
    surf(X2{1},X2{2},F2,'FaceColor','r','EdgeColor','k')
    % surf(X1{1},X1{2},F1,'FaceColor',C(mod(n-1,size(C,1))+1,:),'EdgeColor','k')
    % surf(X2{1},X2{2},F2,'FaceColor',C(mod(n,size(C,1))+1,:),'EdgeColor','k')
    view(-14,12)
    drawnow;
    clf;
end


% figure()
% hold on
% plot(GRID.gblock.x(:,:,1),GRID.gblock.y(:,:,1),'k'); plot(GRID.gblock.x(:,:,1).',GRID.gblock.y(:,:,1).','k'); axis equal
% plot(SUB_GRID.gblock.x(:,:,1),SUB_GRID.gblock.y(:,:,1),'k'); plot(SUB_GRID.gblock.x(:,:,1).',SUB_GRID.gblock.y(:,:,1).','r'); axis equal
% 
% idx_tmp = [3,1,1];
% xp = [CELLS(idx_tmp(1),idx_tmp(2)).fquad.quad_pts];
% vp = [CELLS(idx_tmp(1),idx_tmp(2)).fvec.v];
% quiver(xp(1,:),xp(2,:),vp(1,:),vp(2,:))


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

function [test_fun,test_fun_grad,test_fun_hess] = generate_test_function()

syms x y

% test_fun_sym(x,y) = 1.0 + y + x + y.^2 + x.*y + x.^2;
% test_fun_sym(x,y) = -4.0 + 2.1*y + 0.02*x + 0.006*y.^2 -0.2*x.*y + x.^2;
test_fun_sym(x,y) = (1/2) * ( 999*x.^4 - 888*y.^4 ) - 666*x.*y;
% test_fun_sym(x,y) = (1/12) * ( 999*x.^4 - 888*y.^4 ) ...
%                   + (1/6)  * ( 777*x.^3 .* y + 666*x.*y.^3 ) ...
%                   + (1/4)  * ( 555*x.^2 .* y.^2 ) ...
%                   + (1/6)  * ( 444*x.^3 + 333*y.^3 ) ...
%                   + (1/2)  * ( 222*x.^2 .* y + 111*y.^2 ) ...
%                   + (1/2)  * ( 99*x.^2  + 88 * y.^2 ) ...
%                   + 77*x.*y + 66*x + 55*y;
test_fun = matlabFunction(test_fun_sym);
test_fun_grad = matlabFunction(gradient(test_fun_sym,[x,y]));
test_fun_hess = matlabFunction(hessian(test_fun_sym,[x,y]));

clear x y test_fun_sym

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

function [X,F,delta] = evaluate_taylor_basis_1(T,N,Q,npts)

% get dimension
n_dim = T.n_dim;

X = get_local_interp_grid(Q,n_dim,npts);

[F,delta] = arrayfun(@(x1,y1)T.eval(N,[x1;y1]),X{1},X{2});

end

function [X,F,delta] = evaluate_taylor_basis_2(T,N,Q,npts)

% get dimension
n_dim = T.n_dim;

X = get_local_interp_grid(Q,n_dim,npts);

[F,delta] = arrayfun(@(x1,y1)T.eval(N,[x1;y1]),X{1},X{2});

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