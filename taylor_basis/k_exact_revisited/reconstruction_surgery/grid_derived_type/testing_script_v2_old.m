%% Testing script 2 for the new geometry loading and calculation (07/31/2024)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2022_2024';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
rmpath(genpath(fullfile(parent_dir,'..')))
addpath(genpath(parent_dir));
rmpath(genpath(fullfile(parent_dir,"SENSEI_Translation/")))
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = fullfile(parent_dir,'anisotropic\test_geometry\svf\svf\_svf0017x0009\svf.grd');
grid_t_file = fullfile(parent_dir,'anisotropic\constrained_least_squares\reconstruction_surgery\grid_derived_type\SVF_GRID.mat');
if isfile(grid_t_file)
    load(grid_t_file)
else
    GRID = read_grd_file_to_struct_old( grid_file );
    GRID = grid_type_old(GRID,3);
    save(grid_t_file,"GRID");
end

N = [GRID.gblock(1).imax,GRID.gblock(1).jmax,GRID.gblock(1).kmax];
Nc = max(N - 1,1);
block_info_list(1)                                      = block_info_t_old();
block_info_list(1).block_id                             = 1;
block_info_list(1).Ncells(:)                            = Nc;

block_id   = 1;
% idx        = [24,1,1];
idx        = [9,2,1];
% p0 = polynomial_mod_t_old; p0 = p0.read_from_file('C:\Users\Will\Desktop\test_poly_dump_output.dat');
n_stencil  = 25;
balanced   = true;
print_iter = false;


[ stencil, n_stencil, stencil_cells ] = cell_t_old.build_stencil( block_id, idx, n_stencil, GRID, block_info_list, balanced, print_iter );
stencil_cells = stencil_cells(1:n_stencil);
n_var = 5;
n_dim = 2;
degree = 4;

% soln_stencil = ones(n_stencil,n_var);
inputs = struct();
inputs.gamma        = 1.4;
inputs.radius_inner = 2;
inputs.rho_inner    = 1;
inputs.p_inner      = 1/inputs.gamma;
inputs.mach_inner   = 2.0;

% q_scale = [1.0,75.0,75.0,75.0,100000.0];
q_scale = [1,1,1,1,1];
select = @(M,r) M(r);
qfun1 = @(x)svf_prim(x,inputs);
qfuns = {@(x)select(qfun1(x),1),@(x)select(qfun1(x),2),@(x)select(qfun1(x),3),@(x)select(qfun1(x),4),@(x)select(qfun1(x),5)};

% soln_stencil = calc_soln_in_stencil_quads(n_var,qfuns{var_select},stencil,GRID);
soln_stencil = calc_soln_in_stencil_quads(n_var,qfun1,stencil,GRID);

% soln_stencil = calc_soln_in_stencil_center(n_var,qfun,stencil,GRID);
% set_bc_constraints( this, bound_num, n_var, n_dim, n_constraints, constraint_types, VAR_IDX, constraint_eval )
% GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(1,5,2,1,3000,{[2,3,4]},[0]);
for b = 1:6
    % GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,1,2,1,1000,{1},qfuns(var_select));
    % GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,5,2,1,3000,{[2,3,4]},{@(x)0*x(1)});
    % GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,5,2,5,[1000,1000,1000,1000,1000],{1,2,3,4,5},qfuns);
    GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,4,2,4,[1000,1000,1000,1000],{1,2,3,5},qfuns);
    % GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,1,2,1,[1000],{5},qfuns);
end

% for b = [2,4]
%     GRID.gblock(1) = GRID.gblock(1).set_bc_constraints(b,5,2,1,3000,{[2,3,4]},{@(x)0*x(1)});
% end

neighbor_max_degree = 0;

rec1 = reconstruct_mod_t_old();
rec1 = rec1.setup_interior( n_var, n_dim, degree, stencil, GRID );
rec1.polynomial = rec1.least_squares.solve(soln_stencil,rec1.polynomial,rec1.moments,stencil,GRID);
pfun1 = @(x)rec1.polynomial.evaluate(x);
%
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
plot_var = 2;

Nplot = [25,25,1];

sub_stencil1 = stencil(:,1:9);
sub_stencil2 = sub_stencil1;
% sub_stencil2 = stencil(:,1:1);
[soln_ex,grid_in_stencil] = eval_soln_in_stencil(n_var,qfun1, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil1,GRID);

[soln_ex1,grid_in_stencil1] = eval_soln_in_stencil(n_var,qfun1, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil1,GRID);
[soln_p1,~]               = eval_soln_in_stencil(n_var,pfun1, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil2,GRID);
[soln_p2,~]               = eval_soln_in_stencil(n_var,pfun2, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil2,GRID);
[soln_p3,~]               = eval_soln_in_stencil(n_var,pfun3, q_scale, Nplot(1),Nplot(2),Nplot(3),sub_stencil2,GRID);


k_plane  = 1;

figure();
plot_stencil_cells(plot_var,k_plane,soln_ex1,grid_in_stencil1,'FaceColor','interp','EdgeColor','none');
plot_stencil_cells(plot_var,k_plane,soln_p1,grid_in_stencil1,'FaceColor','r','EdgeColor','none');
plot_stencil_cells(plot_var,k_plane,soln_p2,grid_in_stencil1,'FaceColor','b','EdgeColor','none');
plot_stencil_cells(plot_var,k_plane,soln_p3,grid_in_stencil1,'FaceColor','g','EdgeColor','none');
view(-15,35)

% figure();
% plot_stencil_cells(plot_var,k_plane,soln_ex1,grid_in_stencil1,'FaceColor','interp','EdgeColor','none');
% plot_stencil_cells(plot_var,k_plane,soln_p3,grid_in_stencil1,'FaceColor','g','EdgeColor','none');
% 
% figure();
% plot_stencil_cells(plot_var,k_plane,soln_ex1,grid_in_stencil1,'FaceColor','interp','EdgeColor','none');
% plot_stencil_cells(plot_var,k_plane,soln_p2,grid_in_stencil1,'FaceColor','b','EdgeColor','none');
% 
% figure();
% plot_stencil_cells(plot_var,k_plane,soln_ex1,grid_in_stencil1,'FaceColor','interp','EdgeColor','none');
% plot_stencil_cells(plot_var,k_plane,soln_p1,grid_in_stencil1,'FaceColor','r','EdgeColor','none');
% 
% figure();
% plot_stencil_cells(plot_var,k_plane,soln_ex,grid_in_stencil,'FaceColor','interp','EdgeColor','none');

%%
x = GRID.gblock(1).grid_vars.xi_face_quad(idx(1),idx(2),idx(3)).quad_pts;
n_quad = GRID.gblock(1).grid_vars.xi_face_quad(idx(1),idx(2),idx(3)).n_quad;
normals = GRID.gblock(1).grid_vars.xi_n(idx(1),idx(2),idx(3)).v;
v = zeros(3,n_quad);
vn   = zeros(n_quad,1);
vmag = zeros(n_quad,1);
for n = 1:n_quad
    prim_vals = pfun3(x(:,n));
    v(:,n) = prim_vals(2:4);
    vn(n) = dot( v(:,n),normals(:,n));
    vmag(n) = sqrt( dot(v(:,n),v(:,n)) );
end

soln_stencil_p1 = calc_soln_in_stencil_quads(n_var,pfun1,stencil,GRID);
soln_stencil_p2 = calc_soln_in_stencil_quads(n_var,pfun2,stencil,GRID);
soln_stencil_p3 = calc_soln_in_stencil_quads(n_var,pfun3,stencil,GRID);

m_err1 = soln_stencil_p1(1,:) - soln_stencil(1,:);
m_err2 = soln_stencil_p2(1,:) - soln_stencil(1,:);
m_err3 = soln_stencil_p3(1,:) - soln_stencil(1,:);
%%

fprintf('Error in Conservation of Mean:\n')
fprintf('           density         x-velocity      y-velocity      z-velocity      pressure\n')
num_fmt = '%14.6E';
fprintf(['%s:  ',repmat(sprintf('%s  ',num_fmt),1,n_var-1),sprintf('%s',num_fmt),'\n'],'REC_1',m_err1)
fprintf(['%s:  ',repmat(sprintf('%s  ',num_fmt),1,n_var-1),sprintf('%s',num_fmt),'\n'],'REC_2',m_err2)
fprintf(['%s:  ',repmat(sprintf('%s  ',num_fmt),1,n_var-1),sprintf('%s',num_fmt),'\n'],'REC_3',m_err3)


function plot_stencil_cells(var,k_plane,soln_in_stencil,grid_in_stencil,varargin)
n_stencil = size(soln_in_stencil,1);
hold on;

for n = 1:n_stencil
    X = grid_in_stencil{n,1}(:,:,k_plane);
    Y = grid_in_stencil{n,2}(:,:,k_plane);
    % Z = grid_in_stencil{n,3}(:,:,k_plane);
    V = soln_in_stencil{n,var}(:,:,k_plane);
    % V = V + 1.0e-16*rand(size(V));
    surf(X,Y,V,varargin{:});
    plot_2D_mesh_edge(X,Y,V,'k');
end

X = grid_in_stencil{1,1}(:,:,k_plane);
Y = grid_in_stencil{1,2}(:,:,k_plane);
V = soln_in_stencil{1,var}(:,:,k_plane);
plot_2D_mesh_edge(X,Y,V,'r');

view(2)
axis equal
axis square
axis vis3d
% daspect([1,1,1])
end

function c = contour_stencil_cells(var,k_plane,soln_in_stencil,grid_in_stencil,varargin)
n_stencil = size(soln_in_stencil,1);
hold on;
c = cell(n_stencil+1,1);
for n = 1:n_stencil
    X = grid_in_stencil{n,1}(:,:,k_plane);
    Y = grid_in_stencil{n,2}(:,:,k_plane);
    % Z = grid_in_stencil{n,3}(:,:,k_plane);
    V = soln_in_stencil{n,var}(:,:,k_plane);
    % V = V + 1.0e-12*rand(size(V));
    [~,c{n}] = contourf(X,Y,V,varargin{:});
end
axis equal
X = grid_in_stencil{1,1}(:,:,k_plane);
Y = grid_in_stencil{1,2}(:,:,k_plane);
c{n+1} = plot_2D_mesh_edge(X,Y,0*X,'k');
end


function p = plot_2D_mesh_edge(X,Y,Z,varargin)
hold on;
X1 = X(1:end,1);
X2 = X(end,1:end);
X3 = X(end:-1:1,end);
X4 = X(1,end:-1:end);
Xp = [X1(:);X2(:);X3(:);X4(:)];

Y1 = Y(1:end,1);
Y2 = Y(end,1:end);
Y3 = Y(end:-1:1,end);
Y4 = Y(1,end:-1:end);
Yp = [Y1(:);Y2(:);Y3(:);Y4(:)];

Z1 = Z(1:end,1);
Z2 = Z(end,1:end);
Z3 = Z(end:-1:1,end);
Z4 = Z(1,end:-1:end);
Zp = [Z1(:);Z2(:);Z3(:);Z4(:)];
p = plot3( Xp, Yp, Zp, varargin{:} );
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

function soln_in_stencil = calc_soln_in_stencil_center(n_var,qfun,stencil,grid)
n_stencil = size(stencil,2);
soln_in_stencil = zeros(n_stencil,n_var);
for si = 1:n_stencil
    b   = stencil(1,si);
    idx = stencil(2:4,si);
    soln_in_stencil(si,:) = qfun(grid.gblock(b).grid_vars.cell_c(:,idx(1),idx(2),idx(3)));
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

function q = svf_prim( x, inputs )

    % use set_constants,   only : zero, half, one
    % use fluid_constants, only : gamma, gm1, xg, xgm1
    % ! Non-dimensionalize inputs
    % constructor%radius_inner = constructor%radius_inner * l_ref_grid / l_ref
    % constructor%rho_inner    = constructor%rho_inner    / rho_ref
    % constructor%p_inner      = constructor%p_inner      / ( rho_ref * a_ref**2 )
    % 
    % constructor%a_inner = speed_of_sound( constructor%p_inner,                 &
    %                                       constructor%rho_inner )
    % constructor%u_inner = constructor%mach_inner                               &
    %                     * constructor%rho_inner**( gm1*half )

    q = zeros(5,1);

    radius = sqrt( x(1)^2 + x(2)^2 ); % radius of location

    gamma = inputs.gamma;
    % Density
    q(1) = inputs.rho_inner * ( 1 + (gamma-1)/2*inputs.mach_inner^2 ...
                            * ( 1 - (inputs.radius_inner/radius)^2) ) ...
                            ^(1/(gamma-1));

    % non-dimensionalized ploc = p/(rhoinf*ainf^2)
    ploc = q(1)^gamma/gamma;

    a_inner = sqrt( gamma*inputs.p_inner/inputs.rho_inner);
    u_inner = inputs.mach_inner * inputs.rho_inner^( (gamma-1)/2 );

    % Pressure
    q(5) = ploc*inputs.rho_inner*a_inner^2;

    % x-velocity
    q(2) =  x(2) * a_inner * u_inner * inputs.radius_inner / radius^2;

    % y-velocity
    q(3) = -x(1) * a_inner * u_inner * inputs.radius_inner / radius^2;

end





