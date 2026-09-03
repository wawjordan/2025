%% Experimenting with mapping to reference space (revisited) (08/26/2026)
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
N      = 9;
r      = 4;
o      = 4;
go     = o;
Ng     = r*(N-1)+1;
n_quad = max(2,ceil( 0.5*(o + 1) ));
% Ns     = 2;
Ns     = ceil(1.5*nchoosek(n_dim+o,o));
balanced = true;

[~,~,~,fh1] = my_geomspace_w_refine(Ng,r,0,xmax=1,dx0=0.01/(N-1));
[~,~,~,fh2] = my_geomspace_w_refine(Ng,r,0,xmax=1,dx0=2/(N-1));
% [~,~,~,fh1] = my_geomspace_w_refine(Ng,r,0,xmax=1,dx0=1/(N-1));
% [~,~,~,fh2] = my_geomspace_w_refine(Ng,r,0,xmax=1,dx0=1/(N-1));
x1_map = @(x1,x2,n) x1 + (x2-x1)*fh1(linspace(0,1,n));
x2_map = @(x1,x2,n) x1 + (x2-x1)*fh2(linspace(0,1,n));

% GRID = curv_grid(Ng,Ng,x1_map=x1_map,x2_map=x2_map);
GRID = cart_grid(Ng,Ng,x1_map=x1_map,x2_map=x2_map,theta=pi/3);
FGRID = grid_type(GRID); % fine grid
GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=n_quad,nskip=[r,r,1]);

blk = 1;
idx = [1,1,1];

S = make_stencil(GRID,blk,idx,Ns,balanced);

I = calculate_inertia_tensor( GRID.gblock(S.blk(1)).grid_vars.quad( S.idx(1,1), S.idx(2,1), S.idx(3,1) ) );
[R0,~] = eig(I(1:2,1:2));
R = eye(3);
R(1:2,1:2) = R0;

[Aco,Acon,Mjac] = get_cell_jacobian(FGRID,GRID,blk,idx);
% Acon = Acon/max(abs(Acon),[],'all');

R = get_rotation_matrix(Aco.');

R1 = get_rotation_matrix(Mjac);
R = [0,1,0;1,0,0;0,0,1]*R1;

[n_nodes,xi_vec,x_vec] = get_stencil_centroid_nodes(n_dim,GRID,S);
x0 = x_vec(:,1);
xi0 = xi_vec(:,1);
[n_nodes_eval,xi_vec_eval,x_vec_eval] = get_stencil_quad_nodes(2,n_quad,GRID,S);
A_vec = get_stencil_quad_derivs(2,n_quad,FGRID,GRID,S);
% A_vec = get_stencil_centroid_derivs(FGRID,GRID,S);
% A_vec = A_vec(:,:,1:n_quad^2);
A_vec = A_vec(:,:,1:1);
% n_nodes = n_nodes_eval;
% x_vec = x_vec_eval;
% xi_vec = xi_vec_eval;
% x0 = x_vec(:,1);

n_nodes = n_nodes_eval+1;
x_vec = [x0,x_vec_eval];
xi_vec = [xi0,xi_vec_eval];
x0 = x_vec(:,1);


msk = false(n_nodes,1);
msk(1:1) = true;
% msk(1:n_quad^2) = true;

scale = true;

% fit x(xi)
[x_coefs0,~,cx0] = fit_xi_to_x_map(2,eye(2),xi_vec,x_vec,go,msk,scale);
[x_coefs1,~,cx1] = fit_xi_to_x_map(2,     R,xi_vec,x_vec,go,msk,scale);
[x_coefs2,~,cx2] = fit_xi_to_x_map(2,  Acon,xi_vec,x_vec,go,msk,scale);
x_eval0 = eval_xi_to_x_map(2,go,eye(2),x0,xi_vec_eval,x_coefs0);
x_eval1 = eval_xi_to_x_map(2,go,     R,x0,xi_vec_eval,x_coefs1);
x_eval2 = eval_xi_to_x_map(2,go,  Acon,x0,xi_vec_eval,x_coefs2);

% fit xi(x)
[xi_coefs0,~,c0] = fit_x_to_xi_map(2,eye(2),xi_vec,x_vec,go,msk,scale);
[xi_coefs1,~,c1] = fit_x_to_xi_map(2,     R,xi_vec,x_vec,go,msk,scale);
[xi_coefs2,~,c2] = fit_x_to_xi_map(2,  Acon,xi_vec,x_vec,go,msk,scale);

xi_eval0 = eval_x_to_xi_map(2,go,eye(2),x0,x_vec_eval,xi_coefs0);
xi_eval1 = eval_x_to_xi_map(2,go,     R,x0,x_vec_eval,xi_coefs1);
xi_eval2 = eval_x_to_xi_map(2,go,  Acon,x0,x_vec_eval,xi_coefs2);


% [coefs,M,condM,scale] = fit_x_to_xi_map_w_deriv(2,eye(2),xi_vec,x_vec,Acon*2,go,msk,scale);
% xi_eval2 = eval_x_to_xi_map(2,go,eye(2),x0,x_vec_eval,coefs);

% [coefs,M,condM,scale] = fit_x_to_xi_map_w_deriv(2,Acon,xi_vec,x_vec,Aco.'*8,go,msk,scale);
[coefs,M,condM,scale] = fit_x_to_xi_map_w_deriv(2,Acon,xi_vec,x_vec,A_vec*8,go,msk,scale);
xi_eval3 = eval_x_to_xi_map(2,go,Acon,x0,x_vec_eval,coefs);

N_stencil = S.N;

tiledlayout(1,3)

nexttile
hold on

plot_grid_2D_local(GRID,'k')
plot_grid_2D_local(FGRID,'k:')

for i = 1:N_stencil
    cnt = (i-1)*n_quad^2;
    tmp_x = reshape(x_vec_eval(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    x1 = squeeze(tmp_x(1,:,:));
    x2 = squeeze(tmp_x(2,:,:));
    plot(x1,x2,'k.-')
    plot(x1.',x2.','k.-')

    tmp_x = reshape(x_eval0(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    x1 = squeeze(tmp_x(1,:,:));
    x2 = squeeze(tmp_x(2,:,:));
    plot(x1,x2,'r-')
    plot(x1.',x2.','r-')

    tmp_x = reshape(x_eval1(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    x1 = squeeze(tmp_x(1,:,:));
    x2 = squeeze(tmp_x(2,:,:));
    plot(x1,x2,'b-')
    plot(x1.',x2.','b-')

    tmp_x = reshape(x_eval2(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    x1 = squeeze(tmp_x(1,:,:));
    x2 = squeeze(tmp_x(2,:,:));
    plot(x1,x2,'g-')
    plot(x1.',x2.','g-')
end

% plot(x_eval0(1,:),x_eval0(2,:),'bo')
% plot(x_vec(1,:),x_vec(2,:),'r.')

for i = 2:N_stencil
    % plot_cell_stencil_coords_2D(FGRID,GRID,S,blk,idx,i,r+1,'r-')
    plot_unmapped_stencil_points_2D(FGRID,GRID,S,blk,i,'r-')
end
% plot_cell_stencil_coords_2D(FGRID,GRID,S,blk,idx,1,r+1,'b.-')
plot_unmapped_stencil_points_2D(FGRID,GRID,S,blk,1,'b.-')

% plot_stencil_coords_2D_edges(FGRID,GRID,S,blk,idx,2*r+1,'m','r.-')
% axis equal
% idx = [1,1,1];
% plot_stencil_coords_2D_edges(FGRID,GRID,S,blk,idx,2*r+1,'g','c--')
axis equal
hold off;

nexttile
hold on
N_stencil = S.N;

MM = Acon;
% MM = R;

for i = 2:N_stencil
    plot_mapped_stencil_points_2D(MM,FGRID,GRID,S,blk,idx,i,'r-');
end
plot_mapped_stencil_points_2D(MM,FGRID,GRID,S,blk,idx,1,'b.-')

axis equal
hold off;


nexttile
hold on
for i = 1:N_stencil
    cnt = (i-1)*n_quad^2;
    tmp_xi = reshape(xi_vec_eval(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    xi1 = squeeze(tmp_xi(1,:,:));
    xi2 = squeeze(tmp_xi(2,:,:));
    plot(xi1,xi2,'k.-')
    plot(xi1.',xi2.','k.-')

    % tmp_xi = reshape(xi_eval0(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    % xi1 = squeeze(tmp_xi(1,:,:));
    % xi2 = squeeze(tmp_xi(2,:,:));
    % plot(xi1,xi2,'r.-')
    % plot(xi1.',xi2.','r.-')
    % 
    % tmp_xi = reshape(xi_eval1(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    % xi1 = squeeze(tmp_xi(1,:,:));
    % xi2 = squeeze(tmp_xi(2,:,:));
    % plot(xi1,xi2,'b.-')
    % plot(xi1.',xi2.','b.-')
    % 
    tmp_xi = reshape(xi_eval2(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    xi1 = squeeze(tmp_xi(1,:,:));
    xi2 = squeeze(tmp_xi(2,:,:));
    plot(xi1,xi2,'g.-')
    plot(xi1.',xi2.','g.-')

    tmp_xi = reshape(xi_eval3(:,cnt+1:cnt+n_quad^2),[2,n_quad,n_quad]);
    xi1 = squeeze(tmp_xi(1,:,:));
    xi2 = squeeze(tmp_xi(2,:,:));
    plot(xi1,xi2,'m.-')
    plot(xi1.',xi2.','m.-')
end

axis equal
hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function S = make_all_stencils(GRID,n_stencil,balanced)
% S = struct();
% for blk = 1:GRID.nblocks
%     S.sblock(blk) = struct();
%     IDX = generate_index_array([1 1 1],GRID.gblock(blk).Ncells);
%     S.sblock(blk).stencils = arrayfun( @(idx)make_stencil(GRID,blk,idx{:},n_stencil,balanced), IDX );
% end
% end

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

function [n_nodes,xi_vec,x_vec] = get_stencil_centroid_nodes(n_dim,GRID,S)
n_nodes = S.N;

if ( nargout > 1)
    xi_vec = zeros(n_dim,n_nodes);
    for i = 1:S.N
        off = 2.*(S.idx(:,i)-S.idx(:,1));
        xi_vec(1:n_dim,i) = off(1:n_dim);
    end
end

if ( nargout > 2 )
    x_vec = zeros(n_dim,n_nodes);
    for i = 1:S.N
        x_vec(1:n_dim,i) = GRID.gblock(S.blk(i)).grid_vars.cell_c(1:n_dim,S.idx(1,i),S.idx(2,i),S.idx(3,i));
    end
end
end

function A_vec = get_stencil_centroid_derivs(FGRID,GRID,S)
n_nodes = S.N;
A_vec = zeros(3,3,n_nodes);
for i = 1:S.N
    A_vec(:,:,i) = get_cell_jacobian(FGRID,GRID,S.blk(i),S.idx(:,i).',[0;0;0]).';
end
end

function [n_nodes,xi_vec,x_vec] = get_stencil_quad_nodes(n_dim,n_quad,GRID,S)
switch n_dim
    case(1)
        quad_ref = quad_t.create_quad_ref_1D(n_quad);
    case(2)
        quad_ref = quad_t.create_quad_ref_2D(n_quad);
    case(3)
        quad_ref = quad_t.create_quad_ref_3D(n_quad);
    otherwise
        error('dim must be 1-3')
end
n_nodes = quad_ref.n_quad * S.N;


if ( nargout > 1)
    xi_vec = zeros(n_dim,n_nodes);
    cnt = 0;
    for i = 1:S.N
        off = 2.*(S.idx(:,i)-S.idx(:,1));
        for q = 1:quad_ref.n_quad
            cnt = cnt + 1;
            xi_vec(1:n_dim,cnt) = off(1:n_dim) + quad_ref.quad_pts(1:n_dim,q);
        end
    end
end

if ( nargout > 2 )
    x_vec = zeros(n_dim,n_nodes);
    cnt = 0;
    for i = 1:S.N
        for q = 1:quad_ref.n_quad
            cnt = cnt + 1;
            x_vec(1:n_dim,cnt) = GRID.gblock(S.blk(i)).grid_vars.quad(S.idx(1,i),S.idx(2,i),S.idx(3,i)).quad_pts(1:n_dim,q);
        end
    end
end
end

function A_vec = get_stencil_quad_derivs(n_dim,n_quad,FGRID,GRID,S)
switch n_dim
    case(1)
        quad_ref = quad_t.create_quad_ref_1D(n_quad);
    case(2)
        quad_ref = quad_t.create_quad_ref_2D(n_quad);
    case(3)
        quad_ref = quad_t.create_quad_ref_3D(n_quad);
    otherwise
        error('dim must be 1-3')
end
n_nodes = quad_ref.n_quad * S.N;
A_vec = zeros(3,3,n_nodes);
cnt = 0;
for i = 1:S.N
    for q = 1:quad_ref.n_quad
        cnt = cnt + 1;
        xi_vec = quad_ref.quad_pts(:,q);
        A_vec(:,:,cnt) = get_cell_jacobian(FGRID,GRID,S.blk(i),S.idx(:,i).',xi_vec).';
    end
end
end

function plot_grid_2D_local(GRID,varargin)
for blk = 1:GRID.nblocks
    plot(GRID.gblock(blk).x(:,:,1),GRID.gblock(blk).y(:,:,1),varargin{:});
    plot(GRID.gblock(blk).x(:,:,1).',GRID.gblock(blk).y(:,:,1).',varargin{:});
end
end

% function plot_stencil_coords_2D_edges(FGRID,GRID,S,blk,idx,npts,color,varargin)
% plot_cell_stencil_coords_2D_edges(FGRID,GRID,S,blk,idx,1,npts,varargin{:},'Color',color,'LineStyle','-','Marker','o','MarkerEdgeColor',color)
% for i = 2:S.N
%     plot_cell_stencil_coords_2D_edges(FGRID,GRID,S,blk,idx,i,npts,varargin{:})
% end
% 
% end
% 
% function plot_stencil_coords_2D(FGRID,GRID,S,blk,idx,npts,color,varargin)
% plot_cell_stencil_coords_2D(FGRID,GRID,S,blk,idx,1,npts,varargin{:},'Color',color,'LineStyle','-','Marker','o','MarkerEdgeColor',color)
% for i = 2:S.N
%     plot_cell_stencil_coords_2D(FGRID,GRID,S,blk,idx,i,npts,varargin{:})
% end
% 
% end

% function [M1,M2,M3] = get_cell_jacobian(FGRID,GRID,blk,idx)
% n_skip = GRID.nskip;
% cell_idx = (idx-1).*n_skip+1;
% [xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
% xq = GRID.gblock(blk).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
% xtmp = xtmp - xq(1);
% ytmp = ytmp - xq(2);
% ztmp = ztmp - xq(3);
% L = lagrange_interpolant(size(xtmp,1));
% 
% M1 = L.metrics(xtmp,ytmp,ztmp,[0;0;0]);
% [M2,~] = L.base_vectors(xtmp,ytmp,ztmp,[0;0;0]);
% M2 = 2*M2;
% [xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(GRID.gblock(blk),idx,[1,1,1]);
% xtmp = xtmp - xq(1);
% ytmp = ytmp - xq(2);
% ztmp = ztmp - xq(3);
% M3 = hex_jacobian(xtmp(:),ytmp(:),ztmp(:),0.5,0.5,0.5);
% % M2 = L.base_vectors(xtmp,ytmp,ztmp,xq);
% end
function [Mco,Mcon,Mjac] = get_cell_jacobian(FGRID,GRID,blk,idx,xi)
if (nargin<5)
    xi = [0;0;0];
end
n_skip = GRID.nskip;
cell_idx = (idx-1).*n_skip+1;
[xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
Mjac = hex_jacobian(xtmp(:),ytmp(:),ztmp(:), 0.5, 0.5, 0.5);
xq = GRID.gblock(blk).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
xtmp = xtmp - xq(1);
ytmp = ytmp - xq(2);
ztmp = ztmp - xq(3);
L = lagrange_interpolant(size(xtmp,1));
[Mco,Mcon] = L.base_vectors(xtmp,ytmp,ztmp,xi(:));
Mco  = Mco/2;
Mcon = Mcon.'/2;
end

function plot_cell_stencil_coords_2D(FGRID,GRID,S,blk,idx,i,npts,varargin)

[Xq1,Xq2] = ndgrid(linspace(-1,1,npts),linspace(-1,1,npts));
Xq3 = zeros(npts,npts);

n_skip = GRID.nskip;
cell_idx = (idx-1).*n_skip+1;
[xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
[xfun,yfun,~] = local_grid_interpolant_2D(xtmp,ytmp,ztmp);
off = 2.*(S.idx(:,i)-S.idx(:,1));

X = xfun(Xq1+off(1),Xq2+off(2),Xq3);
Y = yfun(Xq1+off(1),Xq2+off(2),Xq3);
plot(X,Y,varargin{:})
plot(X.',Y.',varargin{:})
end

function plot_unmapped_stencil_points_2D(FGRID,GRID,S,blk,i,varargin)
n_skip = GRID.nskip;
sidx = S.idx(:,i).';
cell_idx = (sidx-1).*n_skip+1;
[xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
X = xtmp(:,:,1);
Y = ytmp(:,:,1);
% Z = ztmp(:,:,1);
plot(X,Y,varargin{:})
plot(X.',Y.',varargin{:})
end

function plot_mapped_stencil_points_2D(M,FGRID,GRID,S,blk,idx,i,varargin)
n_skip = GRID.nskip;
sidx = S.idx(:,i).';
cell_idx = (sidx-1).*n_skip+1;
[xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
xq = GRID.gblock(blk).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
xtmp = xtmp(:,:,1) - xq(1);
ytmp = ytmp(:,:,1) - xq(2);
ztmp = ztmp(:,:,1) - xq(3);
pts = M*[xtmp(:).';ytmp(:).';ztmp(:).'];
X = reshape(pts(1,:),size(xtmp));
Y = reshape(pts(2,:),size(ytmp));
% z = reshape(pts(3,:),size(ztmp));
plot(X,Y,varargin{:})
plot(X.',Y.',varargin{:})
end

% function plot_cell_stencil_coords_2D_edges(FGRID,GRID,S,blk,idx,i,npts,varargin)
% n_skip = GRID.nskip;
% cell_idx = (idx-1).*n_skip+1;
% [xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
% [xfun,yfun,~] = local_grid_interpolant_2D(xtmp,ytmp,ztmp);
% 
% off = 2.*(S.idx(:,i)-S.idx(:,1));
% Xq3 = zeros(npts,1);
% xi_eta = {linspace(-1,1,npts),linspace(-1,1,npts)};
% 
% [Xq1,Xq2] = ndgrid(xi_eta{1}(:),xi_eta{2}(1));
% X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% plot(X,Y,varargin{:})
% 
% [Xq1,Xq2] = ndgrid(xi_eta{1}(:),xi_eta{2}(end));
% X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% plot(X,Y,varargin{:})
% 
% [Xq1,Xq2] = ndgrid(xi_eta{1}(1),xi_eta{2}(:));
% X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% plot(X,Y,varargin{:})
% 
% [Xq1,Xq2] = ndgrid(xi_eta{1}(end),xi_eta{2}(:));
% X = xfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% Y = yfun(Xq1(:)+off(1),Xq2(:)+off(2),Xq3);
% plot(X,Y,varargin{:})
% end

function [xfun,yfun,zfun] = local_grid_interpolant_2D(xtmp,ytmp,ztmp)
L = lagrange_interpolant(size(xtmp,1));
fun = @(xq) L.map_point(xtmp,ytmp,ztmp,xq);
xfun = @(xq1,xq2,xq3) eval_x(fun,xq1,xq2,xq3);
yfun = @(xq1,xq2,xq3) eval_y(fun,xq1,xq2,xq3);
zfun = @(xq1,xq2,xq3) eval_z(fun,xq1,xq2,xq3);
end

% function [xfun,yfun,zfun] = local_grid_interpolant_2D_from_quad(Q,n_quad)
% xtmp = reshape(Q.quad_pts(1,:),n_quad);
% ytmp = reshape(Q.quad_pts(2,:),n_quad);
% ztmp = reshape(Q.quad_pts(3,:),n_quad);
% L = lagrange_interpolant(n_quad(1) + 1);
% fun = @(xq) L.map_point(xtmp,ytmp,ztmp,xq);
% xfun = @(xq1,xq2,xq3) eval_x(fun,xq1,xq2,xq3);
% yfun = @(xq1,xq2,xq3) eval_y(fun,xq1,xq2,xq3);
% zfun = @(xq1,xq2,xq3) eval_z(fun,xq1,xq2,xq3);
% end

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

function GRID = cart_grid(N_x,N_y,varargin)
p = inputParser;
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
validFunctionHandle = @(x) isa(x,'function_handle');
addRequired(p,'N_x',validScalarPosInt);
addRequired(p,'N_y',validScalarPosInt);
addOptional(p,'N_z',nan,validScalarPosInt);
addOptional(p,'x1_map',@(x1,x2,n)linspace(x1,x2,n),validFunctionHandle);
addOptional(p,'x2_map',@(x1,x2,n)linspace(x1,x2,n),validFunctionHandle);
addOptional(p,'x3_map',@(x1,x2,n)linspace(x1,x2,n),validFunctionHandle);
addOptional(p,'theta',0,@(x)isscalar(x));
parse(p,N_x,N_y,varargin{:});

if (nargin(p.Results.x1_map)~=3)
    error('incorrect number of input arguments for x1_map')
end
if (nargin(p.Results.x2_map)~=3)
    error('incorrect number of input arguments for x2_map')
end
if (nargin(p.Results.x3_map)~=3)
    error('incorrect number of input arguments for x3_map')
end
x1_map = p.Results.x1_map;
x2_map = p.Results.x2_map;
x3_map = p.Results.x3_map;
theta  = p.Results.theta;

xi_range   = [0, 1];
eta_range  = [0, 1];
zeta_range = [0, 1];

N_x = p.Results.N_x;
N_y = p.Results.N_y;
N_z = p.Results.N_z;

GRID.dim  = 2;
GRID.imax = N_x;
GRID.jmax = N_y;
xi   = x1_map(xi_range(1),xi_range(2),N_x);
eta  = x2_map(eta_range(1),eta_range(2),N_y);
zeta = 0;
if (~isnan(p.Results.N_z))
    GRID.dim  = 3;
    GRID.kmax = N_z;
    zeta  = x3_map(zeta_range(1),zeta_range(2),N_z);
end

[XI,ETA,ZETA] = ndgrid(xi,eta,zeta);
if (abs(theta)>eps(1))
    R = rotation_matrix_2D(theta);
    [XI,ETA,ZETA] = arrayfun(@(xi,eta,zeta)apply_rotation(xi,eta,zeta,R),XI,ETA,ZETA);
end
GRID.x = XI;
GRID.y = ETA;
if GRID.dim == 3
    GRID.z = ZETA;
end
end

function GRID = curv_grid(N_x,N_y,varargin)
p = inputParser;
validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
validFunctionHandle = @(x) isa(x,'function_handle');
addRequired(p,'N_x',validScalarPosInt);
addRequired(p,'N_y',validScalarPosInt);
addOptional(p,'N_z',nan,validScalarPosInt);
addOptional(p,'x1_map',@(x1,x2,n)linspace(x1,x2,n),validFunctionHandle);
addOptional(p,'x2_map',@(x1,x2,n)linspace(x1,x2,n),validFunctionHandle);
addOptional(p,'x3_map',@(x1,x2,n)linspace(x1,x2,n),validFunctionHandle);
parse(p,N_x,N_y,varargin{:});

if (nargin(p.Results.x1_map)~=3)
    error('incorrect number of input arguments for x1_map')
end
if (nargin(p.Results.x2_map)~=3)
    error('incorrect number of input arguments for x2_map')
end
if (nargin(p.Results.x3_map)~=3)
    error('incorrect number of input arguments for x3_map')
end
x1_map = p.Results.x1_map;
x2_map = p.Results.x2_map;
x3_map = p.Results.x3_map;

xi_range   = [0, 1];
eta_range  = [0, 1];
zeta_range = [0, 1];

N_x = p.Results.N_x;
N_y = p.Results.N_y;
N_z = p.Results.N_z;

GRID.dim  = 2;
GRID.imax = N_x;
GRID.jmax = N_y;
xi   = x1_map(xi_range(1),xi_range(2),N_x);
eta  = x2_map(eta_range(1),eta_range(2),N_y);
zeta = 0;
if (~isnan(p.Results.N_z))
    GRID.dim  = 3;
    GRID.kmax = N_z;
    zeta  = x3_map(zeta_range(1),zeta_range(2),N_z);
end

[XI,ETA,ZETA] = ndgrid(xi,eta,zeta);

GRID.x = xFunc3D(XI,ETA,ZETA);
GRID.y = yFunc3D(XI,ETA,ZETA);

if GRID.dim == 3
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

function vals = eval_poly(coefs,nodevec,degree)
n_dim = size(nodevec,1);
n_nodes = size(nodevec,2);
[exponents,~,~] = get_exponents( n_dim, degree );
n_terms = size(exponents,2);
b = ones(n_terms,1);
vals = zeros(n_nodes,1);
for i = 1:n_nodes
    for j = 1:n_terms
        b(j) = evaluate_monomial( n_dim, exponents,j,nodevec(:,i));
    end
    vals(i) = dot(coefs,b);
end
end

function x_eval = eval_xi_to_x_map(dim,o,A,x0,xi_eval,coefs)
xi_eval = xi_eval(1:dim,:);
n_nodes = size(xi_eval,2);
x0      = x0(1:dim);
A       = A(1:dim,1:dim);
x_eval_tmp = zeros(dim,n_nodes);
for d = 1:dim
    x_eval_tmp(d,:) = eval_poly(coefs(:,d),xi_eval,o);
end
x_eval = A\x_eval_tmp + x0;
end

function xi_eval = eval_x_to_xi_map(dim,o,A,x0,x_eval,coefs)
x_eval = x_eval(1:dim,:);
n_nodes = size(x_eval,2);
x0      = x0(1:dim);
A       = A(1:dim,1:dim);
x_eval = A*(x_eval-x0);
xi_eval = zeros(dim,n_nodes);
for d = 1:dim
    xi_eval(d,:) = eval_poly(coefs(:,d),x_eval,o);
end
end

function [coefs,M,condM,scale] = fit_xi_to_x_map(dim,A,xi_vec,x_vec,o,msk,use_scale)
x0     = x_vec(1:dim,1);
x_vec  = x_vec(1:dim,:) - x0;
xi_vec = xi_vec(1:dim,:);
A      = A(1:dim,1:dim);

M      = computational_transform_matrix(xi_vec,o);
scale = ones(1,size(M,2));
if (use_scale)
    scale = sum( abs(M),1 );
    scale = 1./scale;
end
M = M.*scale;
condM  = cond(M);
x_vec2 = A*x_vec;
% coefs = M\x_vec2.';
n_terms = size(M,2);
coefs  = zeros(n_terms,dim);
for d = 1:dim
    coefs(:,d) = lsqcon_local(M(~msk,:),x_vec2(d,~msk),M(msk,:),x_vec2(d,msk));
    coefs(:,d) = coefs(:,d).*scale.';
end
end

function [coefs,M,condM,scale] = fit_x_to_xi_map(dim,A,xi_vec,x_vec,o,msk,use_scale)
x0     = x_vec(1:dim,1);
x_vec  = x_vec(1:dim,:) - x0;
xi_vec = xi_vec(1:dim,:);
A      = A(1:dim,1:dim);

x_vec2 = A*x_vec;
M      = computational_transform_matrix(x_vec2,o);
scale = ones(1,size(M,2));
if (use_scale)
    scale = sum( abs(M),1 );
    scale = 1./scale;
end
M = M.*scale;
condM  = cond(M);
% coefs = M\xi_vec.';
n_terms = size(M,2);
coefs  = zeros(n_terms,dim);
for d = 1:dim
    coefs(:,d) = lsqcon_local(M(~msk,:),xi_vec(d,~msk),M(msk,:),xi_vec(d,msk));
    coefs(:,d) = coefs(:,d).*scale.';
end
end

function [coefs,M,condM,scale] = fit_x_to_xi_map_w_deriv(dim,A,xi_vec,x_vec,A_vec,o,msk,use_scale)
x0     = x_vec(1:dim,1);
x_vec  = x_vec(1:dim,:) - x0;
xi_vec = xi_vec(1:dim,:);
A      = A(1:dim,1:dim);
A_vec  = A_vec(1:dim,1:dim,:);

for i = 1:size(A_vec,3)
    A_vec(:,:,i) = A*(A_vec(:,:,i).');
end

n_grad = size(A_vec,3);
n_diff_rows = n_grad*dim;

x_vec2 = A*x_vec;
% Md = gradient_constraint_rows(x_vec2(:,1),o);
Md = gradient_constraint_rows(x_vec2(:,1:n_grad),o);
M1 = computational_transform_matrix(x_vec2,o);
M = [Md;M1];
scale = ones(1,size(M,2));
if (use_scale)
    scale = sum( abs(M),1 );
    scale = 1./scale;
end
M = M.*scale;
condM  = cond(M);
Md = M(1:n_diff_rows,:);
M1 = M(n_diff_rows+1:end,:);
% coefs = M\xi_vec.';
n_terms = size(M1,2);
coefs  = zeros(n_terms,dim);
Mcon = M1(msk,:);
Mlsq = [Md;M1(~msk,:)];
for d = 1:dim
    xcon = xi_vec(d,msk);
    tmp = A_vec(:,d,:);
    xlsq = [tmp(:);xi_vec(d,~msk).'];
    coefs(:,d) = lsqcon_local(Mlsq,xlsq,Mcon,xcon);
    coefs(:,d) = coefs(:,d).*scale.';
end
end



function M = computational_transform_matrix(nodevec,degree)
n_dim = size(nodevec,1);
n_nodes = size(nodevec,2);
[exponents,~,~] = get_exponents( n_dim, degree );
n_terms = size(exponents,2);
M = ones(n_nodes,n_terms);
for j = 1:n_terms
    for i = 1:n_nodes
        M(i,j) = evaluate_monomial( n_dim, exponents,j,nodevec(:,i));
    end
end
end

function M = gradient_constraint_rows(nodevec,degree)
n_dim = size(nodevec,1);
n_nodes = size(nodevec,2);
[exponents,~,diff_idx] = get_exponents( n_dim, degree );
n_terms = size(exponents,2);
M = ones(n_nodes*n_dim,n_terms);
for j = 1:n_terms
    cnt = 0;
    for k = 1:n_nodes
        for i = 1:n_dim
            cnt = cnt + 1;
            order = exponents(:,diff_idx(i,1));
            [dval,dcoef,~] = evaluate_monomial_derivative( n_dim, exponents,j,nodevec,order);
            M(cnt,j) = dval*dcoef;
        end
    end
end
end


function [val,coef] = evaluate_monomial( n_dim, exponents, term, x )
val = 1.0;
coef = 1;
for d = 1:n_dim
    for i = exponents(d,term):-1:1
        val = val*x(d);
        coef = coef * i;
    end
end
end

function [dval,dcoef,coef] = evaluate_monomial_derivative( n_dim,     ...
                                                           exponents, ...
                                                           term,      ...
                                                           x,         ...
                                                           order )
dval  = 0.0;
dcoef = 1;
coef  = 1;
if ( any( exponents(:,term) - order(:) < 0 ) ); return; end
dval = 1.0;
for d = 1:n_dim
    for i = exponents(d,term):-1:exponents(d,term)-order(d)+1
        dcoef = dcoef * i;
    end
    for i = exponents(d,term)-order(d):-1:1
        coef = coef * i;
        dval = dval*x(d);
    end
end
end

function [exponents,idx,diff_idx] = get_exponents( n_dim, degree )
n_terms   = nchoosek( n_dim + degree, degree );
exponents = zeros(n_dim,n_terms);
idx       = zeros(degree+1,1);
diff_idx  = zeros(n_dim,n_terms);
cnt = 0;
for curr_total_degree = 0:degree
    nsub(1:n_dim) = curr_total_degree + 1;
    N_full_terms = (curr_total_degree+1)^n_dim;
    for j = 0:N_full_terms
        tmp_exp = global2local(j+1,nsub)-1;
        if ( sum(tmp_exp) == curr_total_degree )
            cnt = cnt + 1;
            exponents(:,cnt) = tmp_exp;
        end
    end
    idx(curr_total_degree+1) = cnt;
end
diff_idx(1:n_dim,1:n_terms) = -1;
if ( degree == 0 ); return; end
for j = 1:idx(degree)
    tmp_exp = exponents(:,j);
    curr_total_degree = sum(tmp_exp);
    cnt = 0;
    for i = idx(curr_total_degree+1)+1:idx(curr_total_degree+2)
        if ( sum( abs(exponents(:,i) - tmp_exp) ) == 1 )
            cnt = cnt +1;
            diff_idx(cnt,j) = i;
        end
    end
end
end

function iSub = global2local(iG,nSub)
nDims = numel(nSub);
iSub = zeros(1,nDims);
if (nDims==1)
    iSub(1) = iG;
    return
end
p = prod(nSub);
iGtmp = iG;
for i = nDims:-1:1
    p = fix( p/nSub(i) );
    iTmp = mod(iGtmp-1,p)+1;
    iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
    iGtmp = iTmp;
end
end

% function iG = local2global(iSub,nSub)
% iSub = iSub(:);
% nSub = nSub(:);
% nDims = numel(iSub);
% p = 1;
% iG = 1;
% for i = 1:nDims
%     iG = iG + ( iSub(i) - 1 )*p;
%     p = p*nSub(i);
% end
% end

function x = lsqcon_local(A, b, B, d)
% compute the QR factorization APb=Q[R,f;0,g]
[Q,R,P] = qr(A);

% n1 = size(A,2);
f = Q \ b(:);

% solve RP^Ty = f  for y
y = R*P.'\f(:);

% solve (P)(R^T)(K^T) = B^T  for K^T
KT = P*R.' \ B.';

% compute the minimum norm solution of Ku = d -  By
u = KT.' \ (d(:) - B*y);

% solve (R)(P^T)z = u  for z
z = R*P.' \ u;

% set x = y + z
x = y + z;

end

function R = rotation_matrix_2D(theta)
R = [cos(theta),-sin(theta),0;...
     sin(theta), cos(theta),0;...
     0,0,1];
end

function [X,Y,Z] = apply_rotation(x,y,z,R)
    tmp = R*[x;y;z];
    X = tmp(1);
    Y = tmp(2);
    Z = tmp(3);
end

function R = get_rotation_matrix(M)
R = zeros(3,3);
ct = M(1,1)/sqrt( M(1,1)^2 + M(2,1)^2 );
st = M(2,1)/sqrt( M(1,1)^2 + M(2,1)^2 );

l1  = vmag(M(:,1));
cp = M(3,1)/l1;
sp = sqrt( M(1,1)^2 + M(2,1)^2 )/l1;

x1 = get_jac_unit_vec(M,1);
x2 = get_jac_unit_vec(M,2);
crx12 = cross( x1, x2 );
cosphi12 = dot( x1, x2 );
sinphi12 = vmag( crx12 );
r2 = (x2 - cosphi12*x1)/sinphi12;
r3 = crx12/sinphi12;
aperp = [-st; +ct; 0];
cb = dot(r2,aperp);
sb = dot(r3,aperp);

R(1,1) = ct*sp;
R(2,1) = st*sp;
R(3,1) = cp;
R(1,2) = -st*cb + ct*cp*sb;
R(2,2) = ct*cb + st*cp*sb;
R(3,2) = -sp*sb;
R(1,3) = -st*sb -ct*cp*cb;
R(2,3) = ct*sb - st*cp*cb;
R(3,3) = sp*cb;

function val = vmag(v)
val = sqrt( sum(v.^2) );
end

function l = get_jac_length(M,dir)
l = vmag(M(:,dir));
end

function vec = get_jac_unit_vec(M,dir)
vec = M(:,dir)/get_jac_length(M,dir);
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
function xc = calculate_centroid( quad )
volume = quad.integrate(ones(1,quad.n_quad));
xc     = quad.integrate(quad.quad_pts)/volume;
end
end