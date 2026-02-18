%% Principal Axis Determination from a quadrature object (01/12/2026)
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

n_dim  = 2;
n_quad = 3;
degree = 3;
alpha = deg2rad(30);
beta  = deg2rad(0);
gamma = deg2rad(0);
sz    = [1,2,3];

v1 = get_hypercube_vertices(sz(1:n_dim));

R = generate_rotation_matrix_3D(alpha,beta,gamma);
R = R(1:n_dim,1:n_dim);

v2 = pagemtimes(R,v1.').';

diff_verts = v1 - pagemldivide(R,v2.').';

GRID1 = create_grid_from_verts(v1,n_quad);
GRID2 = create_grid_from_verts(v2,n_quad);
Q1 = GRID1.gblock.grid_vars.quad;
Q2 = GRID2.gblock.grid_vars.quad;
% basis1 = zero_mean_basis5(n_dim,degree,get_max_cell_extents(v1.',n_dim),Q1);
% basis2 = zero_mean_basis5(n_dim,degree,get_max_cell_extents(v2.',n_dim),Q2);

basis1 = zero_mean_basis5(n_dim,degree,ones(1,n_dim),Q1);
basis2 = zero_mean_basis5(n_dim,degree,ones(1,n_dim),Q2);


% [xp2,yp2,zp2] = L.map_grid(GRID2.gblock.x,GRID2.gblock.y,GRID2.gblock.z,21,21,1);

hold on;
plot_grid_3D(GRID1,n_dim==2,false,'k');
plot_grid_3D(GRID2,n_dim==2,false,'r');
axis equal
plot_basis(GRID1,1,basis1,5,[21,21,1],'EdgeColor','none')
plot_basis(GRID2,1,basis2,5,[21,21,1],'EdgeColor','none')

diff_q_pts = Q1.quad_pts(1:n_dim,:).' - pagemldivide(R,Q2.quad_pts(1:n_dim,:)).';

I1 = calculate_inertia_tensor( Q1 );
I2 = calculate_inertia_tensor( Q2 );

% these are done in the scaled coordinate system...
I1_ = calculate_inertia_tensor_alt( basis1 );
I2_ = calculate_inertia_tensor_alt( basis2 );

diff_eig = sort(eig(I1(1:n_dim,1:n_dim))) - sort(eig(I2(1:n_dim,1:n_dim)));

diff_eig = sort(eig(I1_(1:n_dim,1:n_dim))) - sort(eig(I2_(1:n_dim,1:n_dim)));


[R2,~] = eig(I2(1:n_dim,1:n_dim));
v3 = pagemldivide(R2,v2.').';

% [R2,U,V] = poldecomp(I2);

disp(I1);

function plot_basis(GRID,blk,basis,n,Np,varargin)
[xp,yp,zp] = GRID.gblock(blk).interp.map_grid(GRID.gblock(blk).x,GRID.gblock(blk).y,GRID.gblock(blk).z,Np(1),Np(2),Np(3));
vals = arrayfun(@(x,y,z)basis.eval(n,[x,y,z]),xp,yp,zp);
surf(xp,yp,vals,varargin{:})
end



function h_ref = get_max_cell_extents(node_list,n_dim)
    h_min = min(node_list(1:n_dim,:),[],2);
    h_max = max(node_list(1:n_dim,:),[],2);
    h_ref = 0.5 * (h_max(:) - h_min(:));
end

function GRID = create_grid_from_verts(verts,n_quad)
% verts should be [n_verts,n_dim]
GRID = struct();

n_dim = size(verts,2);
if n_dim>3
    error('does not support arbitrary dimension, you dipshit')
end
shp = repmat(2,1,n_dim);
GRID.x = reshape(verts(:,1),shp);
if n_dim>1
    GRID.y = reshape(verts(:,2),shp);
end
if n_dim>2
    GRID.z = reshape(verts(:,3),shp);
end
GRID = grid_type(GRID,calc_quads=true,nquad=n_quad,ndim=n_dim);
end

function xc = calculate_centroid( quad )
volume = quad.integrate(ones(1,quad.n_quad));
xc     = quad.integrate(quad.quad_pts)/volume;
end

function I = calculate_inertia_tensor_alt( basis )
if basis.degree<2
    error("basis degree needs to be >=2 for use in computing inertia tensor")
end
I = zeros(3,3);
if basis.n_dim == 2
    I(1,1) = basis.moments(6);
    I(2,2) = basis.moments(4);
    I(3,3) = I(1,1) + I(2,2);
    I(1,2) = -basis.moments(5);
    I(2,1) = I(1,2);
elseif basis.n_dim == 3
    I(1,1) = basis.moments(7) + basis.moments(10);
    I(2,2) = basis.moments(5) + basis.moments(10);
    I(3,3) = basis.moments(5) + basis.moments(7);
    I(1,2) = -basis.moments(6);
    I(2,1) = I(1,2);
    I(1,3) = -basis.moments(8);
    I(3,1) = I(1,3);
    I(2,3) = -basis.moments(9);
    I(3,2) = I(2,3);
end
end

function I = calculate_inertia_tensor( quad )
I = zeros(3,3);
xc = calculate_centroid( quad );
x = quad.quad_pts(1,:) - xc(1);
y = quad.quad_pts(2,:) - xc(2);
z = quad.quad_pts(3,:) - xc(3);

I(1,1) = quad.integrate( y.^2 + z.^2 );
I(2,1) = -quad.integrate( x.*y );
I(3,1) = -quad.integrate( x.*z );
I(1,2) = I(2,1);
I(2,2) = quad.integrate( x.^2 + z.^2 );
I(3,2) = -quad.integrate( y.*z );
I(1,3) = I(3,1);
I(2,3) = I(3,2);
I(3,3) = quad.integrate( x.^2 + y.^2 );
end

function [theta, R] = principal_axes(this,quad)
    M_xx = xbar(3, 1);
    M_xy = xbar(2, 2);
    M_yy = xbar(1, 3);
    M = [M_xx, M_xy; M_xy, M_yy];
    [V, D] = eig(M);
    eig_vals = diag(D);
    [~, sort_idx] = sort(eig_vals, 'descend');
    v1 = V(:, sort_idx(1));
    theta = atan2(v1(2), v1(1));
    c = cos(theta); s = sin(theta);
    R = [c, s; -s, c];
end

function [R,U,V] = poldecomp(F)
%POLDECOMP  Performs the polar decomposition of a regular square matrix.
%   [R U V] = POLDECOMP(F) factorizes a non-singular square matrix F such
%   that F=R*U and F=V*R, where
%   U and V are symmetric, positive definite matrices and
%   R is a rotational matrix
%
%   See also EIG, DIAG, REPMAT
% This kind of decomposition is often used in continuum mechanics so it is
% convenient to comment the code that way. From now, we use the matrix 
% formalism of tensors. C is the right Cauchy-Green deformation tensor, 
% F is the deformation tensor, lambda is the stretch.
% Check input
[m,n] = size(F);
if m ~= n
    error('Matrix must be square.');
end
C = F'*F;
[Q0, lambdasquare] = eig(C);
lambda = sqrt(diag((lambdasquare))); % extract the components
% Uinv is the inverse of U and is constructed with the help of Q0. Uinv is
% produced in the same base as F not in the base of its eigenvectors.
Uinv = repmat(1./lambda',size(F,1),1).*Q0*Q0';
% Using the definition, R, U and V can now be calculated
R = F*Uinv;
U = R'*F;
V = F*R';
end