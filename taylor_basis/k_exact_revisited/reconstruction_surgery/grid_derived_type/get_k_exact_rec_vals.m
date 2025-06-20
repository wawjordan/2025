function [xplt,p1,p2] = get_k_exact_rec_vals(x,degree,test_fun,Npts)

n_var = 1;
n_dim = 1;

qfun1 = @(x) test_fun(x(1));
qfuns = {{qfun1}};
q_scale = 1;
Npts_eval = [Npts,1,1];

GRID = struct();
GRID.twod    = false;
GRID.nblocks = 1;
GRID.gblock  = struct();
GRID.gblock.imax = numel(x);
GRID.gblock.jmax = 2;
GRID.gblock.kmax = 2;
[GRID.gblock.x,GRID.gblock.y] = ndgrid( x, [0,1] );
GRID = grid_type_old(GRID,degree);

N = [GRID.gblock(1).imax,GRID.gblock(1).jmax,GRID.gblock(1).kmax];
Nc = max(N - 1,1);
block_info_list(1)                                      = block_info_t_old();
block_info_list(1).block_id                             = 1;
block_info_list(1).Ncells(:)                            = Nc;

% fix quadrature for 1D
for i = 1:Nc(1)+1
    GRID.gblock(1).grid_vars.xi_n(i).n = 1;
    GRID.gblock(1).grid_vars.xi_n(i).v = [1;0;0];
    GRID.gblock(1).grid_vars.xi_face_quad(i).n_quad = 1;
    pts_tmp = GRID.gblock(1).grid_vars.xi_face_quad(i).quad_pts;
    GRID.gblock(1).grid_vars.xi_face_quad(i).quad_pts = [pts_tmp(1,1);0.5;0.5];
    GRID.gblock(1).grid_vars.xi_face_quad(i).quad_wts = 1;
end

% get cell averages of test function
avgs = zeros(Nc(1),1);
for i = 1:Nc(1)
    x_start = GRID.gblock(1).x(i,1,1);
    x_end   = GRID.gblock(1).x(i+1,1,1);
    volume  = x_end - x_start;
    avgs(i) = integral(@(x)test_fun(x),x_start,x_end,'RelTol',1e-16)/volume;
end




for b = 1:2
    GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,1,n_dim,1,1000,{1},qfuns{1});
end


% stencil_info
block_id   = 1;
n_stencil = degree+1+mod(degree,2);
balanced   = true;
print_iter = false;
neighbor_max_degree = 0;

xplt = cell(Nc(1),1);
p1   = cell(Nc(1),1);
p2   = cell(Nc(1),1);

for i = 1:Nc(1)
    idx        = [i,1,1];
    % make stencil
    [ stencil, n_stencil, stencil_cells ] = cell_t_old.build_stencil( block_id, idx, n_stencil, GRID, block_info_list, balanced, print_iter );
    stencil_cells = stencil_cells(1:n_stencil);

    % get cell averages in cells
    % soln_stencil = calc_soln_in_stencil_quads(n_var,qfun1,stencil,GRID);
    idx_tmp = [stencil_cells(:).idx];
    soln_stencil = avgs(idx_tmp(2,:));

    % unconstrained reconstruction
    rec1 = reconstruct_mod_t_old();
    rec1 = rec1.setup_interior( n_var, n_dim, degree, stencil, GRID );
    rec1.polynomial = rec1.least_squares.solve(soln_stencil,rec1.polynomial,rec1.moments,stencil,GRID);
    pfun1 = @(x)rec1.polynomial.evaluate(x);

    rec2 = reconstruct_mod_t_old();
    rec2 = rec2.setup_boundary2( n_var, n_dim, degree, stencil, GRID, neighbor_max_degree, stencil_cells );
    rec2.polynomial = rec2.least_squares.solve(soln_stencil,rec2.polynomial,rec2.moments,stencil,GRID);
    pfun2 = @(x)rec2.polynomial.evaluate(x);
    
    
    
    sub_stencil = stencil(:,1);
    [p1(i),g] = eval_soln_in_stencil(n_var,pfun1, q_scale, Npts_eval(1),Npts_eval(2),Npts_eval(3),sub_stencil,GRID);
    [p2(i),~] = eval_soln_in_stencil(n_var,pfun2, q_scale, Npts_eval(1),Npts_eval(2),Npts_eval(3),sub_stencil,GRID);
    xplt{i} = g{1,1};
end

end


function [soln_in_stencil,grid_in_stencil] = eval_soln_in_stencil(n_var,qfun,scale,N1,N2,N3,stencil,grid)
n_stencil = size(stencil,2);
soln_in_stencil = cell(n_stencil,n_var);
grid_in_stencil = cell(n_stencil,3);
for si = 1:n_stencil
    b   = stencil(1,si);
    idx = stencil(2:4,si);
    idx = idx(:).';
    sz  = [1,1,1];
    [X1,X2,X3] = grid.gblock(b).copy_gblock_nodes(grid.gblock(b),idx,sz);
    [X,Y,Z] = grid.gblock(b).interp.map_grid(X1,X2,X3,N1,N2,N3);
    grid_in_stencil{si,1} = X;
    grid_in_stencil{si,2} = Y;
    grid_in_stencil{si,3} = Z;
    for vi = 1:n_var
        soln_in_stencil{si,vi} = zeros(N1,N2,N3);
    end
    for k = 1:N3
        for j = 1:N2
            for i = 1:N1
                yq = [grid_in_stencil{si,1}(i,j,k),...
                      grid_in_stencil{si,2}(i,j,k),...
                      grid_in_stencil{si,3}(i,j,k)];
                tmp_soln = qfun(yq);
                for vi = 1:n_var
                    soln_in_stencil{si,vi}(i,j,k) = scale(vi)*tmp_soln(vi);
                end
            end
        end
    end
end
end