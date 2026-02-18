%% Experimenting with mapping to reference space (02/17/2026)
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
n_dim    = 2;
N      = 9;
r = 4;
Ng     = r*(N-1)+1;
n_quad = ceil( 0.5*(r + 1) );
% Ns     = 5;
Ns     = ceil(1.5*nchoosek(n_dim+r,r));
balanced = true;

GRID = curv_grid(Ng,Ng);
% GRID = cart_grid(Ng,Ng); GRID.x = GRID.x/r; GRID.y = GRID.y/r;
GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=n_quad,nskip=[r,r,1]);
S    = make_all_stencils(GRID,Ns,balanced);


hold on
plot_grid_2D_local(GRID,'k')
plot_stencil_coords_2D_edges(GRID,S,1,[1,5,1],5,'r','m')
plot_stencil_coords_2D_edges(GRID,S,1,[2,1,1],5,'g','c--')

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = make_all_stencils(GRID,n_stencil,balanced)
S = struct();
for blk = 1:GRID.nblocks
    S.sblock(blk) = struct();
    IDX = generate_index_array([1 1 1],GRID.gblock(blk).Ncells);
    S.sblock(blk).stencils = arrayfun( @(idx)make_stencil(GRID,blk,idx{:},n_stencil,balanced), IDX );
end
end

function stencil = make_stencil(GRID,blk,idx,n_stencil,balanced)
    stencil = struct();
    block_info_list(1)           = block_info_t();
    block_info_list(1).block_id  = blk;
    block_info_list(1).Ncells(:) = GRID.gblock(blk).Ncells;
    print_iter = false;
    [ stencil_idx, n_stencil, ~ ] = cell_t.build_stencil( blk, idx, n_stencil, GRID, block_info_list, balanced, print_iter );
    stencil.blk = stencil_idx(1,:);
    stencil.idx = stencil_idx(2:4,:);
    stencil.N   = n_stencil;
end

function plot_grid_2D_local(GRID,varargin)
for blk = 1:GRID.nblocks
    plot(GRID.gblock(blk).x(:,:,1),GRID.gblock(blk).y(:,:,1),varargin{:});
    plot(GRID.gblock(blk).x(:,:,1).',GRID.gblock(blk).y(:,:,1).',varargin{:});
end
end

function plot_stencil_coords_2D_edges(GRID,S,blk,idx,npts,color,varargin)
plot_cell_stencil_coords_2D_edges(GRID,S,blk,idx,1,npts,varargin{:},'Color',color,'LineStyle','-','Marker','o','MarkerEdgeColor',color)
for i = 2:S.sblock(blk).stencils(idx(1),idx(2),idx(3)).N
    plot_cell_stencil_coords_2D_edges(GRID,S,blk,idx,i,npts,varargin{:})
end

end

function plot_stencil_coords_2D(GRID,S,blk,idx,npts,color,varargin)
plot_cell_stencil_coords_2D(GRID,S,blk,idx,1,npts,varargin{:},'Color',color,'LineStyle','-','Marker','o','MarkerEdgeColor',color)
for i = 2:S.sblock(blk).stencils(idx(1),idx(2),idx(3)).N
    plot_cell_stencil_coords_2D(GRID,S,blk,idx,i,npts,varargin{:})
end

end

function plot_cell_stencil_coords_2D(GRID,S,blk,idx,i,npts,varargin)
n_dim = GRID.gblock(blk).dim;
Q = GRID.gblock(blk).grid_vars.quad(idx(1),idx(2),idx(3));
n_quad = repmat( round( Q.n_quad.^(1/n_dim) ), 1, n_dim );

quad_ref = quad_t.create_quad_ref_2D( n_quad(1) );

lo = min(quad_ref.quad_pts,[],2); lo = 1./lo;
hi = max(quad_ref.quad_pts,[],2); hi = 1./hi;
[Xq1,Xq2] = ndgrid(linspace(lo(1),hi(1),npts),linspace(lo(2),hi(2),npts));
Xq3 = zeros(npts,npts);

[xfun,yfun,~] = local_grid_interpolant_2D(Q,n_quad);

stencil = S.sblock(blk).stencils(idx(1),idx(2),idx(3));
off = (hi-lo).*(stencil.idx(:,i)-stencil.idx(:,1));

X = xfun(Xq1+off(1),Xq2+off(2),Xq3);
Y = yfun(Xq1+off(1),Xq2+off(2),Xq3);
plot(X,Y,varargin{:})
plot(X.',Y.',varargin{:})
end

function plot_cell_stencil_coords_2D_edges(GRID,S,blk,idx,i,npts,varargin)
n_dim = GRID.gblock(blk).dim;
Q = GRID.gblock(blk).grid_vars.quad(idx(1),idx(2),idx(3));
n_quad = repmat( round( Q.n_quad.^(1/n_dim) ), 1, n_dim );

quad_ref = quad_t.create_quad_ref_2D( n_quad(1) );

lo = min(quad_ref.quad_pts,[],2); lo = 1./lo;
hi = max(quad_ref.quad_pts,[],2); hi = 1./hi;

[xfun,yfun,~] = local_grid_interpolant_2D(Q,n_quad);

stencil = S.sblock(blk).stencils(idx(1),idx(2),idx(3));
off = (hi-lo).*(stencil.idx(:,i)-stencil.idx(:,1));
Xq3 = zeros(npts,1);
xi_eta = {linspace(lo(1),hi(1),npts),linspace(lo(2),hi(2),npts)};

[Xq1,Xq2] = ndgrid(xi_eta{1}(:),xi_eta{2}(1));
X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
plot(X,Y,varargin{:})

[Xq1,Xq2] = ndgrid(xi_eta{1}(:),xi_eta{2}(end));
X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
plot(X,Y,varargin{:})

[Xq1,Xq2] = ndgrid(xi_eta{1}(1),xi_eta{2}(:));
X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
plot(X,Y,varargin{:})

[Xq1,Xq2] = ndgrid(xi_eta{1}(end),xi_eta{2}(:));
X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
plot(X,Y,varargin{:})
end

function [xfun,yfun,zfun] = local_grid_interpolant_2D(Q,n_quad)
xtmp = reshape(Q.quad_pts(1,:),n_quad);
ytmp = reshape(Q.quad_pts(2,:),n_quad);
ztmp = reshape(Q.quad_pts(3,:),n_quad);
L = lagrange_interpolant(n_quad(1) + 1);
fun = @(xq) L.map_point(xtmp,ytmp,ztmp,xq);
xfun = @(xq1,xq2,xq3) eval_x(fun,xq1,xq2,xq3);
yfun = @(xq1,xq2,xq3) eval_y(fun,xq1,xq2,xq3);
zfun = @(xq1,xq2,xq3) eval_z(fun,xq1,xq2,xq3);
end

function x = eval_x(fun,xq1,xq2,xq3)
x = arrayfun(@(xq1,xq2,xq3)eval_x_(fun,xq1,xq2,xq3),xq1,xq2,xq3);
end

function y = eval_y(fun,xq1,xq2,xq3)
y = arrayfun(@(xq1,xq2,xq3)eval_y_(fun,xq1,xq2,xq3),xq1,xq2,xq3);
end

function z = eval_z(fun,xq1,xq2,xq3)
z = arrayfun(@(xq1,xq2,xq3)eval_z_(fun,xq1,xq2,xq3),xq1,xq2,xq3);
end

function x = eval_x_(fun,xq1,xq2,xq3)
tmp = fun([xq1,xq2,xq3]);
x = tmp(1);
end

function y = eval_y_(fun,xq1,xq2,xq3)
tmp = fun([xq1,xq2,xq3]);
y = tmp(2);
end

function z = eval_z_(fun,xq1,xq2,xq3)
tmp = fun([xq1,xq2,xq3]);
z = tmp(3);
end

function GRID = cart_grid(N_x,N_y,N_z)
GRID.imax = N_x;
GRID.jmax = N_y;
GRID.dim  = 2;
   
xi  = 0:N_x-1;
eta = 0:N_y-1;

if nargin == 3
    GRID.kmax = N_z;
    GRID.dim = 3;
    zeta = 0:N_z-1;
    [GRID.x,GRID.y,GRID.z] = ndgrid(xi,eta,zeta);
else
    [GRID.x,GRID.y] = ndgrid(xi,eta);
end

end

function GRID = curv_grid(N_x,N_y,N_z)
GRID.imax = N_x;
GRID.jmax = N_y;
x_range = [0, 1];
y_range = [0, 1];
GRID.dim  = 2;
r_val = 1.3;
   
xi  = linspace(x_range(1),x_range(2),N_x);
eta = linspace(y_range(1),y_range(2),N_y);
zeta = 0;
xi  = my_geomspace(N_x,x_range(1),xmax=x_range(2),r=r_val);
eta  = my_geomspace(N_y,y_range(1),xmax=y_range(2),r=r_val);

if nargin == 3
    GRID.kmax = N_z;
    GRID.dim = 3;
    z_range = [0, 1];
    zeta = linspace(z_range(1),z_range(2),N_z);
    zeta  = my_geomspace(N_z,z_range(1),xmax=z_range(2),r=r_val);
end

[XI,ETA,ZETA] = ndgrid(xi,eta,zeta);

GRID.x = xFunc3D(XI,ETA,ZETA);
GRID.y = yFunc3D(XI,ETA,ZETA);

if nargin == 3
    GRID.z = zFunc3D(XI,ETA,ZETA);
end

end

function x = xFunc3D(xi,eta,zeta)
x =  0.10*sin( (2-xi)*pi .*  eta ) ...
  +  0.05*sin(   2   *pi .* zeta ) ...
  +  0.10*eta + 0.05*zeta + 1.00*xi;
end

function y = yFunc3D(xi,eta,zeta)
y = -0.10*sin( (2+eta)*pi .* xi   ) ...
  +  0.05*sin(   2    *pi .* zeta ) ...
  -  0.10*xi + 0.05*zeta + 1.00*eta;
end

function z = zFunc3D(xi,eta,zeta)
z = -0.03*sin(  2    *pi .* xi .* eta ) ...
  -  0.10*xi + 0.05*eta + 1.00*zeta;
end