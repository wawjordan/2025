%% Grid Metric Prototyping v3 (06/29/2026)
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

R1 = 1; R2 = 2; Nr = 17; Nt = 129; Nz = 1;
[x,y,~] = make_annulus_grid(R1,R2,Nr,Nt,Nz);

FD = o_grid_metrics_fd(x,y);
AN = analytic_metrics(x,y);

contourf(x,y,FD.x_xi-AN.x_xi,100,'linecolor','none')
% contourf(x,y,FD.y_xi-AN.y_xi,100,'linecolor','none')
% contourf(x,y,FD.x_eta-AN.x_eta,100,'linecolor','none')
% contourf(x,y,FD.y_eta-AN.y_eta,100,'linecolor','none')
% contourf(x,y,FD.jac_det-AN.jac_det,100,'linecolor','none')
% contourf(x,y,FD.x_xi2-AN.x_xi2,100,'linecolor','none')
% contourf(x,y,FD.x_eta2-AN.x_eta2,100,'linecolor','none')
% contourf(x,y,FD.x_xieta-AN.x_xieta,100,'linecolor','none')
colorbar;
hold on;
plot(x,y,'k')
plot(x.',y.','k')
axis equal

function [X,Y,Z] = make_annulus_grid(R1,R2,Nr,Nt,Nz)
% [R,T,Z] = ndgrid( linspace(R1,R2,Nr), linspace(0,2*pi,Nt), linspace(0,1,Nz));
[T,R,Z] = ndgrid( -linspace(0,2*pi,Nt), linspace(R1,R2,Nr), linspace(0,1,Nz));
X = R.*cos(T);
Y = R.*sin(T);
end

function OUT = analytic_metrics(x,y)
[Nt,Nr] = size(x);
delta_r = 1/(Nr-1);
delta_theta = (-2*pi)/(Nt-1);
OUT = struct();
OUT.x = x;
OUT.y = y;
OUT.r = sqrt(x.^2+y.^2);
OUT.theta = atan2(y,x);
[dxdtheta, dydtheta, dxdr, dydr, ...
 dxdtheta2, dydtheta2, dxdr2, dydr2, ...
 dxdrdtheta, dydrdtheta] = polar_jacobian(OUT.r,OUT.theta);

OUT.jac_det = -OUT.r * (delta_r*delta_theta);
OUT.x_xi   = dxdtheta    * delta_theta;
OUT.y_xi   = dydtheta    * delta_theta;
OUT.x_eta  = dxdr        * delta_r;
OUT.y_eta  = dydr        * delta_r;
OUT.x_xi2  = dxdtheta2   * delta_theta^2;
OUT.y_xi2  = dydtheta2   * delta_theta^2;
OUT.x_eta2  = dxdr2      * delta_r^2;
OUT.y_eta2  = dydr2      * delta_r^2;
OUT.x_xieta = dxdrdtheta * (delta_r * delta_theta);
OUT.y_xieta = dydrdtheta * (delta_r * delta_theta);
end

function [dxdtheta, dydtheta, dxdr, dydr, ...
          dxdtheta2, dydtheta2, dxdr2, dydr2, ...
          dxdrdtheta, dydrdtheta ] = polar_jacobian(r,theta,~)
dxdr     =        cos(theta);
dydr     =        sin(theta);
dxdtheta =    -r.*sin(theta);
dydtheta =     r.*cos(theta);
dxdr2    = zeros(size(theta));
dydr2    = zeros(size(theta));
dxdtheta2 =   -r.*cos(theta);
dydtheta2 =   -r.*sin(theta);
dxdrdtheta =     -sin(theta);
dydrdtheta =      cos(theta);
end

function OUT = o_grid_metrics_fd(x,y)
OUT = struct();
% [Ni,Nj] = size(x);
OUT.x = x;
OUT.y = y;
[OUT.x_xi,OUT.y_xi,OUT.x_eta,OUT.y_eta] = metrics_2D(x,y);
[OUT.x_xi2,OUT.y_xi2,OUT.x_eta2,OUT.y_eta2,OUT.x_xieta,OUT.y_xieta] = second_order_metrics_2D(x,y);
[~,~,~,OUT.jac_det] = cov_g_2D(OUT.x_xi,OUT.y_xi,OUT.x_eta,OUT.y_eta);
end

function [x_xi,y_xi,x_eta,y_eta] = metrics_2D(x,y)
x_xi = central2_diff_periodic1(x,1);
y_xi = central2_diff_periodic1(y,1);
x_eta = central2_diff_non_periodic1(x,2);
y_eta = central2_diff_non_periodic1(y,2);
end

function [x_xi2,y_xi2,x_eta2,y_eta2,x_xieta,y_xieta] = second_order_metrics_2D(x,y)
x_xi2 = central2_diff_periodic2(x,1);
y_xi2 = central2_diff_periodic2(y,1);
x_eta2 = central2_diff_non_periodic2(x,2);
y_eta2 = central2_diff_non_periodic2(y,2);
x_xieta = central_cross_diff(x,1);
y_xieta = central_cross_diff(y,1);
end

function [g11,g22,g12,J] = cov_g_2D(x_xi,y_xi,x_eta,y_eta)
% covariant metric tensor for 2D curvilinear space
g11 = x_xi.^2 + y_xi.^2;
g22 = x_eta.^2 + y_eta.^2;
g12 = x_xi.*x_eta + y_xi.*y_eta;
J   = x_xi.*y_eta - x_eta.*y_xi;
end

function [g11,g22,g12,J] = con_g_2D(x_xi,y_xi,x_eta,y_eta)
% contravariant metric tensor for 2D curvilinear space
% vectorized, space-saving form
J = x_xi.*y_eta - x_eta.*y_xi;
g0inv = 1./(J.^2);
g11 =  (x_eta.^2 + y_eta.^2).*g0inv;
g12 = -(x_xi.*x_eta + y_xi.*y_eta).*g0inv;
g22 = (x_xi.^2 + y_xi.^2).*g0inv;
end

% function [g11,g22,g12] = cov2con_g_2D(g_11,g_22,g_12,J)
% % contravariant metric tensor for 2D curvilinear space
% % computed from covariant metric tensor
% g0inv = 1./(J.^2);
% g11 =  g_22.*g0inv;
% g12 = -g_12.*g0inv;
% g22 =  g_11.*g0inv;
% end

function dA = central2_diff_periodic1(A,dim)
[NI,NJ] = size(A);
dA = 0*A;
if dim == 1
    for j = 1:NJ
        for i = 2:NI-1
            dA(i,j) = 0.5*( A(i+1,j) - A(i-1,j) );
        end
    end
    for j = 1:NJ
        dA(NI,j) = 0.5*( A(2,j) - A(NI-1,j) );
        dA(1,j) = dA(NI,j);
    end
elseif dim == 2
    for j = 2:NJ-1
        for i = 1:NI
            dA(i,j) = 0.5*( A(i,j+1) - A(i,j-1) );
        end
    end
    for i = 1:NI
        dA(i,NJ) = 0.5*( A(i,2) - A(i,NJ-1) );
        dA(i,1) = dA(i,NJ);
    end
else
    error('dim must be 1 or 2')
end
end
function dA = central2_diff_non_periodic1(A,dim)
[NI,NJ] = size(A);
dA = 0*A;
if dim == 1
    for j = 1:NJ
        for i = 2:NI-1
            dA(i,j) = 0.5*( A(i+1,j) - A(i-1,j) );
        end
    end
    for j = 1:NJ
        dA(NI,j) = 0.5*( 3*A(NI,j) - 4*A(NI-1,j) + A(NI-2,j));
        dA(1,j) = -0.5*( 3*A(1,j) - 4*A(2,j) + A(3,j));
    end
elseif dim == 2
    for j = 2:NJ-1
        for i = 1:NI
            dA(i,j) = 0.5*( A(i,j+1) - A(i,j-1) );
        end
    end
    for i = 1:NI
        dA(i,NJ) = 0.5*( 3*A(i,NJ) - 4*A(i,NJ-1) + A(i,NJ-2));
        dA(i,1) = -0.5*( 3*A(i,1) - 4*A(i,2) + A(i,3));
    end
else
    error('dim must be 1 or 2')
end
end
function dA = central2_diff_periodic2(A,dim)
dA = central2_diff_periodic1(A,dim);
dA = central2_diff_periodic1(dA,dim);
end
function dA = central2_diff_non_periodic2(A,dim)
dA = central2_diff_non_periodic1(A,dim);
dA = central2_diff_non_periodic1(dA,dim);
end

function dA = central12_diff_cross(A,dim)
% dim is dimension of periodicity
dim1 = dim;
dim2 = 1 + (dim1==1);
dA = central2_diff_periodic1(A,dim1);
dA = central2_diff_non_periodic1(dA,dim2);
end

function dA = central21_diff_cross(A,dim)
% dim is dimension of periodicity
dim1 = dim;
dim2 = 1 + (dim1==1);
dA = central2_diff_non_periodic1(A,dim2);
dA = central2_diff_periodic1(dA,dim1);
end

function dA = central_cross_diff(A,dim)
dA = 0.5*(central12_diff_cross(A,dim)+central21_diff_cross(A,dim));
end