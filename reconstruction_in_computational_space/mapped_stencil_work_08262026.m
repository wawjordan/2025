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
r      = 1;
o      = 2;
Ng     = r*(N-1)+1;
n_quad = ceil( 0.5*(o + 1) );
% Ns     = 2;
Ns     = ceil(1.5*nchoosek(n_dim+o,o));
balanced = true;

[~,~,~,fh1] = my_geomspace_w_refine(Ng,r,0,xmax=1,dx0=1.0/(N-1));
[~,~,~,fh2] = my_geomspace_w_refine(Ng,r,0,xmax=1,dx0=1.0/(N-1));
x1_map = @(x1,x2,n) x1 + (x2-x1)*fh1(linspace(0,1,n));
x2_map = @(x1,x2,n) x1 + (x2-x1)*fh2(linspace(0,1,n));

GRID = curv_grid(Ng,Ng,x1_map=x1_map,x2_map=x2_map);
% GRID = cart_grid(Ng,Ng); GRID.x = GRID.x/r; GRID.y = GRID.y/r;
FGRID = grid_type(GRID); % fine grid
GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=n_quad,nskip=[r,r,1]);

blk = 1;
idx = [6,6,1];

S = make_stencil(GRID,blk,idx,Ns,balanced);

[n_nodes,xi_vec,x_vec] = get_stencil_quad_nodes(2,n_quad,GRID,S);

M = computational_transform_matrix(xi_vec,o);

coefs = M\x_vec.';

x_vals = eval_poly(coefs(:,1),xi_vec,o);
y_vals = eval_poly(coefs(:,2),xi_vec,o);



N_stencil = S.N;

tiledlayout(1,2)

nexttile
hold on
plot_grid_2D_local(GRID,'k')
plot_grid_2D_local(FGRID,'k:')
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
A = get_cell_jacobian(FGRID,GRID,blk,idx);
N_stencil = S.N;
for i = 2:N_stencil
    plot_mapped_stencil_points_2D(A,FGRID,GRID,S,blk,idx,i,'r-');
end
plot_mapped_stencil_points_2D(A,FGRID,GRID,S,blk,idx,1,'b.-')

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
function M = get_cell_jacobian(FGRID,GRID,blk,idx)
n_skip = GRID.nskip;
cell_idx = (idx-1).*n_skip+1;
[xtmp,ytmp,ztmp] = GRID.gblock(blk).copy_gblock_nodes(FGRID.gblock(blk),cell_idx,n_skip);
xq = GRID.gblock(blk).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
xtmp = xtmp - xq(1);
ytmp = ytmp - xq(2);
ztmp = ztmp - xq(3);
L = lagrange_interpolant(size(xtmp,1));
[~,M] = L.base_vectors(xtmp,ytmp,ztmp,[0;0;0]);
M = M.'/2;
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