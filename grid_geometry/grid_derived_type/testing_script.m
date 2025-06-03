%% Testing script for the new geometry loading and calculation (04/01/2025)
clc; clear; close all;

agglom=false;
load_file=true;

grid_file = 'C:\Users\Will\Desktop\kt.grd';
if agglom
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_agglom.mat';
    % grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID_agglom_cart.mat';
else
    grid_t_file = 'C:\Users\Will\Documents\MATLAB\VT_Research\2025\grid_geometry\grid_derived_type\GRID.mat';
end
if isfile(grid_t_file)&&load_file
    load(grid_t_file)
else
    GRID = read_grd_file_to_struct( grid_file );
    if agglom
        GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=2,nskip=[2,2,1]);
    else
        GRID = grid_type(GRID,calc_quads=true,nquad=2);
    end
    save(grid_t_file,"GRID");
end

blk  = 1;
idx  = [65,2,1];
dim  = 2;

face_idx = 3;


[center,nodes,face_quads,normals,nbor_center] = get_cell_geom(GRID,blk,idx,dim);

[x_gamma,x0,~] = get_optimal_face_quads(face_quads{face_idx},normals{face_idx});

[M1,beta] = simple_correction_matrix(GRID,blk,idx,dim);
D1_1 = one_exact_derivative(GRID,blk,idx,dim,[1,1,1,1,1]);
D1_1 = one_exact_derivative_op(GRID,blk,idx,dim);

hold on
plot_cell_nodes(nodes,2,'k');
plot_cell_nodes(face_quads{face_idx}.quad_pts,2,'k.');
plot(x_gamma(1),x_gamma(2),'ro');
plot(x0(1),x0(2),'go');

plot(center(1),center(2),'rx')
for i = 1:length(nbor_center)
    plot(nbor_center{i}(1),nbor_center{i}(2),'kx')
end



% p = plot_2D_mesh_edge( GRID.gblock(1).x, GRID.gblock(1).y, GRID.gblock(1).z );

% func = @(x1,x2) arrayfun(@(x1,x2)face_point_func(face_quads{face_idx},normals{face_idx},[x1,x2,0]),x1,x2);
% fsurf(func,[-0.05 0.05 -0.05 0.05])
% function [x_gamma] = get_optimal_face_quads(face_quad,normals)
% 
% x_gamma = zeros(3,1);
% S0 = face_quad.integrate(normals.v);
% 
% % x_1 = face_quad.integrate( face_quad.quad_pts.*normals.v(1,:) );
% % x_2 = face_quad.integrate( face_quad.quad_pts.*normals.v(2,:) );
% % x_3 = face_quad.integrate( face_quad.quad_pts.*normals.v(3,:) );
% 
% x_gamma(1) = ( S0(1) * dot( face_quad.quad_wts, face_quad.quad_pts(1,:).*normals.v(1,:) ) ...
%              + S0(2) * dot( face_quad.quad_wts, face_quad.quad_pts(1,:).*normals.v(2,:) ) ...
%              + S0(3) * dot( face_quad.quad_wts, face_quad.quad_pts(1,:).*normals.v(3,:) ) )...
%              / ( S0(1)^2 + S0(2)^2 + S0(3)^2 );
% x_gamma(2) = ( S0(1) * dot( face_quad.quad_wts, face_quad.quad_pts(2,:).*normals.v(1,:) ) ...
%              + S0(2) * dot( face_quad.quad_wts, face_quad.quad_pts(2,:).*normals.v(2,:) ) ...
%              + S0(3) * dot( face_quad.quad_wts, face_quad.quad_pts(2,:).*normals.v(3,:) ) )...
%              / ( S0(1)^2 + S0(2)^2 + S0(3)^2 );
% x_gamma(3) = ( S0(1) * dot( face_quad.quad_wts, face_quad.quad_pts(3,:).*normals.v(1,:) ) ...
%              + S0(2) * dot( face_quad.quad_wts, face_quad.quad_pts(3,:).*normals.v(2,:) ) ...
%              + S0(3) * dot( face_quad.quad_wts, face_quad.quad_pts(3,:).*normals.v(3,:) ) )...
%              / ( S0(1)^2 + S0(2)^2 + S0(3)^2 );
% end

% function func = face_point_func(face_quad,normals,x)
% 
% S0 = face_quad.integrate(normals.v);
% 
% S1 = @(x) [ face_quad.integrate( (face_quad.quad_pts(1,:) - x(1)).*normals.v(1,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(1,:) - x(1)).*normals.v(2,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(1,:) - x(1)).*normals.v(3,:) ); ...
%             face_quad.integrate( (face_quad.quad_pts(2,:) - x(2)).*normals.v(1,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(2,:) - x(2)).*normals.v(2,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(2,:) - x(2)).*normals.v(3,:) ); ...
%             face_quad.integrate( (face_quad.quad_pts(3,:) - x(3)).*normals.v(1,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(3,:) - x(3)).*normals.v(2,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(3,:) - x(3)).*normals.v(3,:) ) ];
% func = sum( (S1(x)*S0).^2 );
% 
% 
% end
% 
% function [x_gamma,x0] = get_optimal_face_quads_2(face_quad,normals)
% 
% S0 = face_quad.integrate(normals.v);
% 
% S1 = @(x) [ face_quad.integrate( (face_quad.quad_pts(1,:) - x(1)).*normals.v(1,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(1,:) - x(1)).*normals.v(2,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(1,:) - x(1)).*normals.v(3,:) ); ...
%             face_quad.integrate( (face_quad.quad_pts(2,:) - x(2)).*normals.v(1,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(2,:) - x(2)).*normals.v(2,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(2,:) - x(2)).*normals.v(3,:) ); ...
%             face_quad.integrate( (face_quad.quad_pts(3,:) - x(3)).*normals.v(1,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(3,:) - x(3)).*normals.v(2,:) ), ...
%             face_quad.integrate( (face_quad.quad_pts(3,:) - x(3)).*normals.v(3,:) ) ];
% obj_fun = @(x) sum( (S1(x)*S0).^2 );
% 
% x0 = face_quad.integrate( (face_quad.quad_pts )./ face_quad.integrate(ones(1,face_quad.n_quad)));
% 
% 
% options = optimset('Display','Iter','TolFun',1e-5,'TolX',1e-5);
% x_gamma = fminsearch(obj_fun,x0,options);
% end

function [x_gamma,x0,S0] = get_optimal_face_quads(face_quad,normals)

x_gamma = zeros(3,1);
S0 = face_quad.integrate(normals.v);

for i = 1:3
    x_gamma(i) = dot( S0, face_quad.integrate( face_quad.quad_pts(i,:).*normals.v ) );
end

x_gamma = x_gamma / dot(S0,S0);

x0 = face_quad.integrate( (face_quad.quad_pts )./ face_quad.integrate(ones(1,face_quad.n_quad)));

end

function beta = get_weight_coefficient(x_gamma,x_j,x_k)
beta = norm(x_j - x_gamma,2) / ( norm(x_j - x_gamma,2) + norm(x_k - x_gamma,2) );
end

function [M1,beta] = simple_correction_matrix(grid,blk,idx,dim)
    M1 = zeros(dim,dim);
    [center,~,face_quads,normals,nbor_center] = get_cell_geom(grid,blk,idx,dim);
    n = length(nbor_center);
    beta = zeros(1,n);
    AIJ  =  cell(1,n);
    for k = 1:n
        [x_gamma,~,AIJ{k}] = get_optimal_face_quads(face_quads{k},normals{k});
        beta(k) = get_weight_coefficient(x_gamma,center,nbor_center{k});
        % M1 = M1 + beta(k) * ( nbor_center{k}(1:dim) - center(1:dim) ) * AIJ{k}(1:dim).';
        M1 = M1 + beta(k) * AIJ{k}(1:dim) * ( nbor_center{k}(1:dim) - center(1:dim) ).';
    end
    M1 = M1/grid.gblock(blk).grid_vars.volume(idx(1),idx(2),idx(3));
end

function D1_1 = one_exact_derivative_op(grid,blk,idx,dim)
    [M1,beta] = simple_correction_matrix(grid,blk,idx,dim);
    [~,~,face_quads,normals,~] = get_cell_geom(grid,blk,idx,dim);
    n = length(normals);
    D1 = zeros(dim,n+1);
    for k = 1:n
        AJK = face_quads{k}.integrate(normals{k}.v);
        D1(:,1)   = D1(:,1) + (1 - beta(k))* AJK(1:dim);
        D1(:,1+k) = beta(k) * AJK(1:dim);
    end
    D1 = D1 / grid.gblock(blk).grid_vars.volume(idx(1),idx(2),idx(3));
    D1_1 = M1 \ D1;
end


function D1_1 = one_exact_derivative(grid,blk,idx,dim,avgs)
    [M1,beta] = simple_correction_matrix(grid,blk,idx,dim);
    [~,~,face_quads,normals,~] = get_cell_geom(grid,blk,idx,dim);
    n = length(normals);
    D1 = zeros(dim,1);

    for k = 1:n
        AJK = face_quads{k}.integrate(normals{k}.v);
        D1 = D1 + ( beta(k)*avgs(1+k) + (1 - beta(k))*avgs(1) ) * AJK(1:dim);
    end
    D1 = D1 / grid.gblock(blk).grid_vars.volume(idx(1),idx(2),idx(3));
    D1_1 = M1 \ D1;
end

function [center,nodes,face_quads,normals,nbor_center] = get_cell_geom(GRID,blk,idx,dim)

bnd_min = [1;1;1];
bnd_max = [ GRID.gblock(blk).imax - 1;...
            GRID.gblock(blk).jmax - 1;...
            GRID.gblock(blk).kmax - 1 ];

center    = GRID.gblock(blk).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
node_var={'x','y','z'};
quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
norm_var={'xi_n','eta_n','zeta_n'};
if dim == 2
    nodes      = zeros(3,4);
    face_quads = cell(1,4);
    normals    = cell(1,4);
    nbor_center = cell(1,4);
    off   = [0,1,1,0; ...
             0,0,1,1; ...
             0,0,0,0];
    for i = 1:3
        for j = 1:4
            nodes(i,j) = GRID.gblock(blk).(node_var{i})( idx(1) + off(1,j),...
                                                         idx(2) + off(2,j),...
                                                         idx(3) + off(3,j) );
        end
    end
    for i = 1:4
        face_idx =  fix((i-1)/2)+1;      % [1,1,2,2]
        offset   =  mod(i+1,2);          % [0,1]
        factor   =  2*offset-1;          % [-1, 1]
        tmp_idx  =  idx;
        tmp_idx(face_idx) = tmp_idx(face_idx)+offset;
        face_quads{i} = GRID.gblock(blk).grid_vars.(quad_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}    = GRID.gblock(blk).grid_vars.(norm_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}.v = normals{i}.v * factor;
        tmp_idx  =  idx;
        tmp_idx(face_idx) = tmp_idx(face_idx)+factor;
        if in_bounds(tmp_idx,bnd_min,bnd_max)
            nbor_center{i}  = GRID.gblock(blk).grid_vars.cell_c(:,tmp_idx(1),tmp_idx(2),tmp_idx(3));
        end
    end
else
    nodes      = zeros(3,8);
    face_quads = cell(1,6);
    normals    = cell(1,6);
    nbor_center = cell(1,6);
    off   = [0,1,1,0,0,1,1,0; ...
             0,0,1,1,0,0,1,1; ...
             0,0,0,0,1,1,1,1];
    for i = 1:3
        for j = 1:8
            nodes(i,j) = GRID.gblock(blk).(node_var{i})( idx(1) + off(1,j),...
                                                         idx(2) + off(2,j),...
                                                         idx(3) + off(3,j) );
        end
    end
    for i = 1:6
        face_idx =  fix((i-1)/2)+1; % [1,1,2,2,3,3]
        offset   =  mod(i+1,2);     % [0,1]
        factor   =  2*offset-1;     % [-1, 1]
        tmp_idx  =  idx;
        tmp_idx(face_idx) = tmp_idx(face_idx)+offset;
        face_quads{i} = GRID.gblock(blk).grid_vars.(quad_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}    = GRID.gblock(blk).grid_vars.(norm_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}.v = normals{i}.v * factor;
        if in_bounds(tmp_idx,bnd_min,bnd_max)
            nbor_center{i}  = GRID.gblock(blk).grid_vars.cell_c(:,tmp_idx(1),tmp_idx(2),tmp_idx(3));
        end
    end
end
end

function val = in_bounds(idx,bnd_min,bnd_max)
val = ( all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:)) || ...
        all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:)) );
end

function p = plot_cell_nodes(nodes,dim,varargin)

if dim == 2
    p = plot( [nodes(1,:),nodes(1,1)], ...
              [nodes(2,:),nodes(2,1)], varargin{:} );
else
    p = plot3( [nodes(1,:),nodes(1,1)], ...
               [nodes(2,:),nodes(2,1)], ...
               [nodes(3,:),nodes(3,1)], varargin{:} );
end
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