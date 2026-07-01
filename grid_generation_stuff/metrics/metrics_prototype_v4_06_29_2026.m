%% Grid Metric Prototyping v4 (06/29/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
grid_file = fullfile(parent_dir,'kt.grd');
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

GRID = read_grd_file_to_struct(grid_file);
x = GRID.gblock.x;
y = GRID.gblock.y;
i1 = 49;
offset_i = 10;
offset_j = 10;
x2 = x(i1-offset_i:end+1-i1+offset_i,1:offset_j);
y2 = y(i1-offset_i:end+1-i1+offset_i,1:offset_j);
i1 = 1+offset_i;
FD = c_grid_jmin_metrics_fd(x2,y2,i1);

contourf(x2,y2,FD.y_eta,100,'linecolor','none')
colorbar;
hold on;
plot(x,y,'k')
plot(x.',y.','k')
axis equal
xlim([-1,2])
ylim([-1.5,1.5])


function OUT = c_grid_jmin_metrics_fd(x,y,i1)
OUT = struct();
% [Ni,Nj] = size(x);
OUT.x = x;
OUT.y = y;
[OUT.x_xi,OUT.y_xi,OUT.x_eta,OUT.y_eta] = metrics_2D_wake_cut_jmin(x,y,i1);
[OUT.x_xi2,OUT.y_xi2,OUT.x_eta2,OUT.y_eta2,OUT.x_xieta,OUT.y_xieta] = second_order_metrics_2D_wake_cut_jmin(x,y,i1);
[~,~,~,OUT.jac_det] = cov_g_2D(OUT.x_xi,OUT.y_xi,OUT.x_eta,OUT.y_eta);
end

function [x_xi,y_xi,x_eta,y_eta] = metrics_2D_wake_cut_jmin(x,y,i1)
x_xi = central2_diff_jmin_wake_cut1(x,1,i1);
y_xi = central2_diff_jmin_wake_cut1(y,1,i1);
x_eta = central2_diff_jmin_wake_cut1(x,2,i1);
y_eta = central2_diff_jmin_wake_cut1(y,2,i1);
end

function [x_xi2,y_xi2,x_eta2,y_eta2,x_xieta,y_xieta] = second_order_metrics_2D_wake_cut_jmin(x,y,i1)
x_xi2 = central2_diff_jmin_wake_cut2(x,1,i1);
y_xi2 = central2_diff_jmin_wake_cut2(y,1,i1);
x_eta2 = central2_diff_jmin_wake_cut2(x,2,i1);
y_eta2 = central2_diff_jmin_wake_cut2(y,2,i1);
x_xieta = central_cross_diff(x,i1);
y_xieta = central_cross_diff(y,i1);
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

% function dA = central2_diff_periodic1(A,dim)
% [NI,NJ] = size(A);
% dA = 0*A;
% if dim == 1
%     for j = 1:NJ
%         for i = 2:NI-1
%             dA(i,j) = 0.5*( A(i+1,j) - A(i-1,j) );
%         end
%     end
%     for j = 1:NJ
%         dA(NI,j) = 0.5*( A(2,j) - A(NI-1,j) );
%         dA(1,j) = dA(NI,j);
%     end
% elseif dim == 2
%     for j = 2:NJ-1
%         for i = 1:NI
%             dA(i,j) = 0.5*( A(i,j+1) - A(i,j-1) );
%         end
%     end
%     for i = 1:NI
%         dA(i,NJ) = 0.5*( A(i,2) - A(i,NJ-1) );
%         dA(i,1) = dA(i,NJ);
%     end
% else
%     error('dim must be 1 or 2')
% end
% end

function dA = central2_diff_jmin_wake_cut1(A,dim,i1)
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
    % for j = 1:NJ
    %     dA(NI,j) = 0.5*( A(2,j) - A(NI-1,j) );
    %     dA(1,j) = dA(NI,j);
    % end
elseif dim == 2
    for j = 2:NJ-1
        for i = 1:NI
            dA(i,j) = 0.5*( A(i,j+1) - A(i,j-1) );
        end
    end
    for i = 1:NI
        dA(i,NJ) = 0.5*( 3*A(i,NJ) - 4*A(i,NJ-1) + A(i,NJ-2));
    end
    for i = 1:i1
        dA(i,1) = 0.5*( A(i,2) - A(NI+1-i,2) );
        dA(NI+1-i) = -dA(i,1);
    end
else
    error('dim must be 1 or 2')
end
end

function dA = central2_diff_jmin_wake_cut2(A,dim,i1)
dA = central2_diff_jmin_wake_cut1(A,dim,i1);
dA = central2_diff_jmin_wake_cut1(dA,dim,i1);
end

function dA = central_diff_cross(A,dim,i1)
% dim is dimension of periodicity
dim1 = dim;
dim2 = 1 + (dim1==1);
dA = central2_diff_jmin_wake_cut1(A,dim1,i1);
dA = central2_diff_jmin_wake_cut1(dA,dim2,i1);
end
function dA = central_cross_diff(A,i1)
dA = 0.5*(central_diff_cross(A,1,i1)+central_diff_cross(A,2,i1));
end