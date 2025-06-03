%% Testing the taylor_shape_functions type (04/30/2025)
clc; clear; close all;

agglom=false;
load_file=true;
cart=true;

grid_file = 'C:\Users\Will\Desktop\kt.grd';

if cart
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_cart.mat';
elseif agglom
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_agglom.mat';
else
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID.mat';
    % grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_quad_5.mat';
end

if isfile(grid_t_file)&&load_file
    load(grid_t_file)
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

blk     = 1;
idx     = [5,5,1];
% idx     = [65,4,1];
% idx     = [33,4,1];
dim     = 2;
degree  = 2;
n_vars  = 1;
bnd_min = ones(dim,1);
bnd_max = GRID.gblock(blk).Ncells(1:dim);


syms x y

% test_fun_sym(x,y) = 1.0 + y + x + y.^2 + x.*y + x.^2;
test_fun_sym(x,y) = -4.0 + 2.1*y + 0.02*x + 0.006*y.^2 -0.2*x.*y + x.^2;
test_fun = matlabFunction(test_fun_sym);
test_fun_grad = matlabFunction(gradient(test_fun_sym,[x,y]));
test_fun_hess = matlabFunction(hessian(test_fun_sym,[x,y]));
% clear x y test_fun_sym
% test_fun = @(x,y) 1.0 + y + x + y.^2 + x.*y + x.^2;
% test_fun = @(x,y) sin(50*x).*cos(50*y);

hold on
idx_low  = [-2,-2,0] + idx;
idx_high = [ 2, 2,0] + idx;

idx_low  = [1;1;1];
idx_high = GRID.gblock(blk).Ncells(:);

SUB_GRID = GRID.subset_grid(1,idx_low,idx_high);

% plot(GRID.gblock.x(:,:,1),GRID.gblock.y(:,:,1),'k'); plot(GRID.gblock.x(:,:,1).',GRID.gblock.y(:,:,1).','k'); axis equal
% plot(SUB_GRID.gblock.x(:,:,1),SUB_GRID.gblock.y(:,:,1),'k'); plot(SUB_GRID.gblock.x(:,:,1).',SUB_GRID.gblock.y(:,:,1).','k'); axis equal

[TS,QS] = get_face_stencil_stencil_shape_functions(GRID,blk,idx,degree);
Cs = get_linear_reconstruction_of_function_over_stencil_stencil(TS,QS,test_fun);
Cs2 = get_exact_linear_rec_of_function_over_stencil_stencil(TS,QS,test_fun,test_fun_grad);

% get cell averages
avgs   = get_cell_avg_of_function(QS(1).Q,test_fun);
dvals2 = get_derivative_approximation(TS(1).T,avgs);

coefs1 = get_quadratic_reconstruction_in_cell(TS(1).T,Cs);
coefs2 = get_quadratic_reconstruction_in_cell(TS(1).T,Cs2);
coefs3 = get_quadratic_reconstruction_in_cell_hard_coded(TS(1).T,QS(1).Q,Cs);
coefs4 = get_quadratic_reconstruction_in_cell_hard_coded(TS(1).T,QS(1).Q,Cs2);

coefs5 = get_exact_quadratic_reconstruction_in_cell(TS(1).T,QS(1).Q,test_fun,test_fun_grad,test_fun_hess);


coefs_plot = coefs5;

figure(2)
hold on;
C = colororder(gcf);
for n = 1:1%length(QS(1).Q)
    [X,F1] = evaluate_function_on_interp_grid(test_fun,QS(1).Q(n),21);
    % surf(X{1},X{2},F1,'FaceColor',C(mod(n-1,size(C,1))+1,:),'EdgeColor','none')
    % plot3([X{1}(:,1);X{1}(end,:).';flip(X{1}(:,end));flip(X{1}(1,:)).'],...
    %       [X{2}(:,1);X{2}(end,:).';flip(X{2}(:,end));flip(X{2}(1,:)).'],...
    %       [   F1(:,1);   F1(end,:).';flip(   F1(:,end));flip(   F1(1,:)).'],'k')
    % 
    % [X,F] = evaluate_reconstruction_on_interp_grid(Cs(n).coefs,TS(1).T(n),QS(1).Q(n),21);
    % surf(X{1},X{2},F,'FaceColor','r','EdgeColor','none')
    % plot3([X{1}(:,1);X{1}(end,:).';flip(X{1}(:,end));flip(X{1}(1,:)).'],...
    %       [X{2}(:,1);X{2}(end,:).';flip(X{2}(:,end));flip(X{2}(1,:)).'],...
    %       [   F(:,1);   F(end,:).';flip(   F(:,end));flip(   F(1,:)).'],'k')
end
[~,F] = evaluate_reconstruction_on_interp_grid(coefs_plot,TS(1).T(1),QS(1).Q(1),21);
surf(X{1},X{2},F-F1,'FaceColor','g','EdgeColor','k')

% [X,F] = evaluate_reconstruction_on_interp_grid(coefs1,TS(1).T(1),QS(1).Q(1),21);
% surf(X{1},X{2},F-0*F1,'FaceColor','b','EdgeColor','none')
% plot3([X{1}(:,1);X{1}(end,:).';flip(X{1}(:,end));flip(X{1}(1,:)).'],...
%       [X{2}(:,1);X{2}(end,:).';flip(X{2}(:,end));flip(X{2}(1,:)).'],...
%       [   F(:,1);   F(end,:).';flip(   F(:,end));flip(   F(1,:)).'],'k')

% M2 = [ Ts(1).build_reconstruction_matrix(Ts(2)); ...
%        Ts(1).build_reconstruction_matrix(Ts(3)); ...
%        Ts(1).build_reconstruction_matrix(Ts(4)) ];
% R2 = [ Ts(1).build_reconstruction_rhs(Ts(2),coefs,nbor_coefs);

view(-14,12)

figure(3)
hold on;
for i = 1:length(TS(1).T)
    C = colororder(gcf);
    for n = 1:TS(1).T(1).n_terms
        [X,F,~] = evaluate_shape_function(TS(1).T(i),n,QS(1).Q(i),21);
        surf(X{1},X{2},F,'FaceColor',C(mod(n-1,size(C,1))+1,:),'EdgeColor','k')
        plot3([X{1}(:,1);X{1}(end,:).';flip(X{1}(:,end));flip(X{1}(1,:)).'],...
              [X{2}(:,1);X{2}(end,:).';flip(X{2}(:,end));flip(X{2}(1,:)).'],...
              [   F(:,1);   F(end,:).';flip(   F(:,end));flip(   F(1,:)).'],'k')
    end
end
view(-14,12)