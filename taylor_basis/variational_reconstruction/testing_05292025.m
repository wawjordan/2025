%% Testing the taylor_shape_functions type (04/30/2025)
clc; clear; close all;

agglom=false;
load_file=true;
cart=true;
GRID = load_gen_grid_for_testing(agglom,load_file,cart);

blk     = 1;
idx     = [5,5,1];
dim     = 2;
degree  = 3;
n_vars  = 1;
[test_fun,~,~] = generate_test_function();

idx_low  = [-2,-2,0] + idx;
idx_high = [ 2, 2,0] + idx;

SUB_GRID = GRID.subset_grid(1,idx_low,idx_high);

CELLS = set_up_cells(SUB_GRID,1,[1,1,1],[5,5,1],degree,n_vars);
CELLS = arrayfun(@(CELLS)CELLS.set_cell_avg({test_fun}),CELLS);
CELLS2 = set_up_cell_var_recs(SUB_GRID,1,[1,1,1],[5,5,1],degree,{test_fun});

[A1,~,~] = CELLS(1).calc_cell_matrices_explicit(CELLS,[1,1,1]);
[A12,~,~] = CELLS(1).calc_cell_matrices(CELLS,[1,1,1]);
A2 = CELLS2(1).get_self_LHS_CELLS(CELLS2,[1,1,1]);
A22 = CELLS2(1).get_self_LHS_CELLS_NEW(CELLS2,[1,1,1]);


A_LHS1 = CELLS(1).construct_full_LHS_matrix(CELLS);
A_LHS2 = var_rec_t.construct_full_LHS_matrix(CELLS2);

RHS1 = CELLS(1).construct_full_RHS(CELLS);
RHS2 = var_rec_t.construct_full_RHS(CELLS2);

coeffs1 = A_LHS1\RHS1;
coeffs2 = A_LHS2\RHS2;


point = CELLS(1).fquad(1).quad_pts(:,1);

T1 = CELLS(1).taylor;
T2 = CELLS2(1).basis;
l = 5;
[tmp1,delta1] = T1.eval_derivative( 5, point, T1.exponents(:,2) );
[tmp2,delta2] = T2.deval( 5, point, T2.terms(2).exponents );

Q1 = CELLS(1).quad;
Q2 = CELLS2(1).quad;

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


figure(2)
hold on
plot(GRID.gblock.x(:,:,1),GRID.gblock.y(:,:,1),'k'); plot(GRID.gblock.x(:,:,1).',GRID.gblock.y(:,:,1).','k'); axis equal
plot(SUB_GRID.gblock.x(:,:,1),SUB_GRID.gblock.y(:,:,1),'k'); plot(SUB_GRID.gblock.x(:,:,1).',SUB_GRID.gblock.y(:,:,1).','r'); axis equal

idx_tmp = [3,1,1];
xp = [CELLS(idx_tmp(1),idx_tmp(2)).fquad.quad_pts];
vp = [CELLS(idx_tmp(1),idx_tmp(2)).fvec.v];
quiver(xp(1,:),xp(2,:),vp(1,:),vp(2,:))

figure(2)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GRID = load_gen_grid_for_testing(agglom,load_file,cart)

grid_file = 'C:\Users\Will\Desktop\kt.grd';

if cart
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_cart.mat';
elseif agglom
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_agglom.mat';
else
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID.mat';
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
        GRID = grid_type(GRID,calc_quads=true,nquad=5);
    end
    save(grid_t_file,"GRID");
end

end

function [test_fun,test_fun_grad,test_fun_hess] = generate_test_function()

syms x y

% test_fun_sym(x,y) = 1.0 + y + x + y.^2 + x.*y + x.^2;
test_fun_sym(x,y) = -4.0 + 2.1*y + 0.02*x + 0.006*y.^2 -0.2*x.*y + x.^2;
test_fun = matlabFunction(test_fun_sym);
test_fun_grad = matlabFunction(gradient(test_fun_sym,[x,y]));
test_fun_hess = matlabFunction(hessian(test_fun_sym,[x,y]));

clear x y test_fun_sym

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