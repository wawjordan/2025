%% Testing script for 1D k-exact (06/20/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_var = 1;
n_dim = 1;
degree = 4;

% make a "1D" grid
x = linspace(-1,1,11);
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

% stencil_info
block_id   = 1;
idx        = [3,1,1];
n_stencil = degree+1+mod(degree,2);
balanced   = true;
print_iter = false;

% make stencil
[ stencil, n_stencil, stencil_cells ] = cell_t_old.build_stencil( block_id, idx, n_stencil, GRID, block_info_list, balanced, print_iter );
stencil_cells = stencil_cells(1:n_stencil);

qfun1 = @(x) sinc(x(1));
% qfuns = {@(x)select(qfun1(x),1)};
qfuns = {{qfun1}};
q_scale = 1;
soln_stencil = calc_soln_in_stencil_quads(n_var,qfun1,stencil,GRID);
for b = 1:2
    GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,1,n_dim,1,1000,{1},qfuns{1});
end
neighbor_max_degree = 0;

rec1 = reconstruct_mod_t_old();
rec1 = rec1.setup_interior( n_var, n_dim, degree, stencil, GRID );
rec1.polynomial = rec1.least_squares.solve(soln_stencil,rec1.polynomial,rec1.moments,stencil,GRID);
pfun1 = @(x)rec1.polynomial.evaluate(x);
mean_factor = 1;
boundary_factor = 1000;
rec2 = reconstruct_mod_t_old();
rec2 = rec2.setup_boundary( n_var, n_dim, degree, stencil, GRID, neighbor_max_degree, stencil_cells, mean_factor, boundary_factor );
rec2.polynomial = rec2.least_squares.solve(soln_stencil,rec2.polynomial,rec2.moments,stencil,GRID);
pfun2 = @(x)rec2.polynomial.evaluate(x);

rec3 = reconstruct_mod_t_old();
rec3 = rec3.setup_boundary2( n_var, n_dim, degree, stencil, GRID, neighbor_max_degree, stencil_cells );
rec3.polynomial = rec3.least_squares.solve(soln_stencil,rec3.polynomial,rec3.moments,stencil,GRID);
pfun3 = @(x)rec3.polynomial.evaluate(x);


%%
plot_var = 1;

Nplot = [21,1,1];

sub_stencil1 = stencil(:,1:end);
sub_stencil2 = sub_stencil1;
% sub_stencil2 = stencil(:,1:1);
[soln_ex,grid_in_stencil] = eval_soln_in_stencil(n_var,qfun1, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil1,GRID);

[soln_ex1,grid_in_stencil1] = eval_soln_in_stencil(n_var,qfun1, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil1,GRID);
[soln_p1,~]               = eval_soln_in_stencil(n_var,pfun1, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil2,GRID);
[soln_p2,~]               = eval_soln_in_stencil(n_var,pfun2, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil2,GRID);
[soln_p3,~]               = eval_soln_in_stencil(n_var,pfun3, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil2,GRID);


k_plane  = 1;

figure();
plot_stencil_cells(plot_var,k_plane,soln_ex1,grid_in_stencil1);
plot_stencil_cells(plot_var,k_plane,soln_p1,grid_in_stencil1,'r');
plot_stencil_cells(plot_var,k_plane,soln_p2,grid_in_stencil1,'b');
plot_stencil_cells(plot_var,k_plane,soln_p3,grid_in_stencil1,'g');


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


function soln_in_stencil = calc_soln_in_stencil_quads(n_var,qfun,stencil,grid)
n_stencil = size(stencil,2);
soln_in_stencil = zeros(n_stencil,n_var);
for si = 1:n_stencil
    b   = stencil(1,si);
    idx = stencil(2:4,si);
    quad = grid.gblock(b).grid_vars.quad(idx(1),idx(2),idx(3));
    vol  = grid.gblock(b).grid_vars.volume(idx(1),idx(2),idx(3));
    tmp_soln = zeros(n_var,quad.n_quad);
    for n = 1:quad.n_quad
        
        tmp_soln(:,n) = qfun(quad.quad_pts(:,n));
    end
    soln_in_stencil(si,:) = quad.integrate(tmp_soln)/vol;
end
end

function plot_stencil_cells(var,k_plane,soln_in_stencil,grid_in_stencil,varargin)
n_stencil = size(soln_in_stencil,1);
hold on;

for n = 1:n_stencil
    X = grid_in_stencil{n,1}(:,:,k_plane);
    % Y = grid_in_stencil{n,2}(:,:,k_plane);
    % Z = grid_in_stencil{n,3}(:,:,k_plane);
    V = soln_in_stencil{n,var}(:,:,k_plane);
    % V = V + 1.0e-16*rand(size(V));
    % surf(X,Y,V,varargin{:});
    plot(X,V,varargin{:});
    % plot_2D_mesh_edge(X,Y,V,'k');
end

% X = grid_in_stencil{1,1}(:,:,k_plane);
% Y = grid_in_stencil{1,2}(:,:,k_plane);
% V = soln_in_stencil{1,var}(:,:,k_plane);
% plot_2D_mesh_edge(X,Y,V,'r');
% 
% view(2)
% axis equal
% axis square
% axis vis3d
% daspect([1,1,1])
end