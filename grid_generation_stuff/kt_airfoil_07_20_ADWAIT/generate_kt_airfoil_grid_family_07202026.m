%% KT airfoil test case (ADWAIT) (07/20/2026)
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
dstr = char(datetime('now',"Format",'uuuu-MM-dd''_''HH-mm-ss'));

levels = [1:6];
n_levels = length(levels);
airfoil_inputs = struct();
airfoil_inputs.epsilon = 0.1;
airfoil_inputs.kappa   = 0.2;
airfoil_inputs.tau     = deg2rad(10.0);
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = -10; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;

n_refine   = 2^4;
n_body_pts = 129;
n_wake_pts = 129;
jmax_in    = 129;
% n_body_pts = 257;
% n_wake_pts = 129;
% jmax_in    = 129;
AR_LE      =  1.0;
AR_TE      =  1.0;
delta_LE   = 0.125;
delta_TE   = 0.125;
spline_order = 5;


% n_body_pts = 65;
% n_wake_pts = 65;
% jmax_in    = 65;
% F1 = transfinite_scalar_interpolant_2D( G1 = 1.0, ...
%                                         G2 = @(t) smooth_ramps(t,[ 0.0, 0.25,  1.000 ], ...
%                                                                  [ 1.0, 0.9999, 0.98 ]), ...
%                                         G3 = 0.98,...
%                                         G4 = @(t) smooth_ramps(t,[ 0.0,  0.25, 1.000 ], ...
%                                                                  [ 1.0, 0.9999, 0.98 ]) );
% F2 = transfinite_scalar_interpolant_2D( G1 = 0.5, ...
%                                         G2 = @(t) smooth_ramps(t,[ 0.00, 0.25,  0.75, 1.0 ], ...
%                                                                  [ 0.50, 2.00, 1.00, 0.5 ]), ...
%                                         G3 = 0.5, ...
%                                         G4 = @(t) smooth_ramps(t,[ 0.00,  0.25, 0.75, 1.0 ], ...
%                                                                  [ 0.50, 2.00, 1.00, 0.5 ]) );

n_body_pts = 129;
n_wake_pts = 129;
jmax_in    = 129;
F1 = transfinite_scalar_interpolant_2D( G1 = 1.0, ...
                                        G2 = @(t) smooth_ramps(t,[ 0.0, 0.25,  1.000 ], ...
                                                                 [ 1.0, 0.9999, 0.98 ]), ...
                                        G3 = 0.98,...
                                        G4 = @(t) smooth_ramps(t,[ 0.0,  0.25, 1.000 ], ...
                                                                 [ 1.0, 0.9999, 0.98 ]) );
F2 = transfinite_scalar_interpolant_2D( G1 = 0.5, ...
                                        G2 = @(t) smooth_ramps(t,[ 0.00, 0.25,  0.75, 1.0 ], ...
                                                                 [ 0.50, 4.00, 4.00, 0.5 ]), ...
                                        G3 = 0.5, ...
                                        G4 = @(t) smooth_ramps(t,[ 0.00,  0.25, 0.75, 1.0 ], ...
                                                                 [ 0.50, 4.00, 4.00, 0.5 ]) );
% n_body_pts = 257;
% n_wake_pts = 257;
% jmax_in    = 257;
% F1 = transfinite_scalar_interpolant_2D( G1 = 1.0, ...
%                                         G2 = @(t) smooth_ramps(t,[ 0.0, 0.25,  1.000 ], ...
%                                                                  [ 1.0, 0.9999, 0.99 ]), ...
%                                         G3 = 0.99,...
%                                         G4 = @(t) smooth_ramps(t,[ 0.0,  0.25, 1.000 ], ...
%                                                                  [ 1.0, 0.9999, 0.99 ]) );
% F2 = transfinite_scalar_interpolant_2D( G1 = 0.5, ...
%                                         G2 = @(t) smooth_ramps(t,[ 0.00, 0.25,  0.75, 1.0 ], ...
%                                                                  [ 0.50, 4.00, 4.00, 0.5 ]), ...
%                                         G3 = 0.5, ...
%                                         G4 = @(t) smooth_ramps(t,[ 0.00,  0.25, 0.75, 1.0 ], ...
%                                                                  [ 0.50, 4.00, 4.00, 0.5 ]) );

imax = n_body_pts + 2*(n_wake_pts-1);
scjmax = @(i,j) F1.eval((i-1)/(imax-1),(j-1)/(jmax_in-1));
alpham = @(i,j) F2.eval((i-1)/(imax-1),(j-1)/(jmax_in-1));
[OUT,airfoil] = generate_kt_airfoil_grid_07_21_2026(airfoil_inputs,n_body_pts,n_wake_pts,jmax_in,AR_LE,AR_TE,delta_LE,delta_TE,spline_order,alpham,scjmax);


% [OUT,airfoil] = generate_kt_airfoil_grid(airfoil_inputs,n_body_pts,n_wake_pts,jmax_in,AR_LE,AR_TE,delta_LE,delta_TE,spline_order);

folder = 'C:\Users\wajordan\Downloads\kt_for_transfer_AR_1_1_TE_10_d-01_ALPHA_5';
out_folder = fullfile(folder,'\grids\');
jobfmt = '_kt%0.4dx%0.4d';
% save(['OUT_STRUCT_',dstr],"OUT");

skip1 = 2;
skip2 = 2;
hold on
plot(OUT.base_grid.x(1:2^skip2:end,1:2^skip1:end),OUT.base_grid.y(1:2^skip2:end,1:2^skip1:end),'k')
plot(OUT.base_grid.x(1:2^skip2:end,1:2^skip1:end).',OUT.base_grid.y(1:2^skip2:end,1:2^skip1:end).','k')
axis equal


figure()
plot_edge_length_ratio(OUT.base_grid.x(1:2^skip2:end,1:2^skip1:end),OUT.base_grid.y(1:2^skip2:end,1:2^skip1:end),1,'jet',1);
xlim([-1,2])
ylim([-1.5,1.5])

plot_edge_length_ratio(OUT.base_grid.x(1:2^skip2:end,1:2^skip1:end),OUT.base_grid.y(1:2^skip2:end,1:2^skip1:end),2,'jet',1);
xlim([-1,2])
ylim([-1.5,1.5])

% L1 = edge_length_ratio_2D(x2,y2,1);
% L2 = edge_length_ratio_2D(x2,y2,2);


imax = n_refine*(OUT.base_grid.imax-1)+1;
jmax = n_refine*(OUT.base_grid.jmax-1)+1;
% imax = 311;
% jmax = 111;
[GRID,i1,i2] = get_fine_grid(OUT,airfoil,imax,jmax);

% skip1 = 0;
% skip2 = 0;
hold on;
plot(GRID.x(1:2^skip2:end,1:2^skip1:end),GRID.y(1:2^skip2:end,1:2^skip1:end),'r');
plot(GRID.x(1:2^skip1:end,1:2^skip2:end).',GRID.y(1:2^skip1:end,1:2^skip2:end).','r');
% plot(GRID.x(1:2^skip2:end,1:2^skip1:10),GRID.y(1:2^skip2:end,1:2^skip1:10),'r');
% plot(GRID.x(1:2^skip1:end,1:2^skip2:10).',GRID.y(1:2^skip1:end,1:2^skip2:10).','r');
axis equal
xlim([-1,2])
ylim([-1.5,1.5])

xlim([0.9868    0.9950])
ylim([-0.0036    0.0046])

bc_id_list  = [-1,-1,201,-200,-200,-200];

Nfine = [GRID.imax,GRID.jmax];
for j = 1:n_levels
    s = 2^(levels(j)-1);
    fprintf('writing level %d of %d\n',j,n_levels)
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [out_folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\kt'];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
    Ni(1) = 1;
    Ni(2) = (i1 - 1)/s + 1;
    Ni(3) = (i2 - 1)/s + 1;
    Ni(4) = (GRID.imax - 1)/s + 1;

    Nj(1) = 1;
    Nj(2) = (GRID.jmax - 1)/s + 1;

    N = [GRID2.imax,GRID2.jmax];

    idx_list{1} =          [ Ni(1), Ni(2), Nj(1), Nj(1) ];
    additional_inputs{1} = [ Ni(4), Ni(3), Nj(1), Nj(1),  1 ];
    idx_list{2} =          [ Ni(3), Ni(4), Nj(1), Nj(1) ];
    additional_inputs{2} = [ Ni(2), Ni(1), Nj(1), Nj(1),  1 ];
    idx_list{3} = [ Ni(2), Ni(3),  Nj(1), Nj(1) ];
    idx_list{4} = [ Ni(1), Ni(4),  Nj(2), Nj(2) ];
    idx_list{5} = [ Ni(1), Ni(1),  Nj(1), Nj(2) ];
    idx_list{6} = [ Ni(4), Ni(4),  Nj(1), Nj(2) ];
    additional_inputs{3} = [];
    additional_inputs{4} = [];
    additional_inputs{5} = [];
    additional_inputs{6} = [];
    SENSEI_BC_write(N,idx_list,bc_id_list,additional_inputs,[filename,'.bc'])
end

function hfig = plot_edge_length_ratio(x,y,dir,cmap,linewidth)
hfig = figure();
set(gca,'Color','k');
colormap(cmap);
L = edge_length_ratio_2D(x,y,dir);
x = padarray(x,[1,1],nan,'both');
y = padarray(y,[1,1],nan,'both');
L = padarray(L,[1,1],nan,'both');
hold on;
patch(x,y,L,'edgecolor','interp','Marker','.','LineWidth',linewidth);
patch(x.',y.',L.','edgecolor','interp','LineWidth',linewidth);
axis equal
colorbar;
hold off
end

function [GRID,i1,i2] = get_fine_grid(OUT,airfoil,imax,jmax)
GRID = struct();
GRID.imax = imax;
GRID.jmax = jmax;
ni = (GRID.imax-1)/(OUT.base_grid.imax-1);
nj = (GRID.jmax-1)/(OUT.base_grid.jmax-1);
i1 = (OUT.base_grid.i1-1)*ni+1;
i2 = (OUT.base_grid.i2-1)*ni+1;
% integer refinement/coarsening
if mod(ni,1)<10*eps(1) && mod(nj,1)<10*eps(1)
    xi_vec  = linspace(0,1,imax);
    eta_vec = linspace(0,1,jmax);
else
    % we need to make sure that there are points on the trailing edge
    i1 = round(i1);
    i2 = round(i2);
    off  = ([5,5]-1)/(imax-1);
    xi1 = [0,i1-1,i2-1,imax-1];
    xi1 = xi1/xi1(end);
    xi2 = [0,(OUT.base_grid.i1-1)/(OUT.base_grid.imax-1),(OUT.base_grid.i2-1)/(OUT.base_grid.imax-1),1.0];
    [fh,~,~] = hermite_blended_piecewise_linear(xi1,xi2,off,off);
    xi_vec  = fh(linspace(0,1,imax));
    eta_vec = linspace(0,1,jmax);
end
GRID.x = fnval(OUT.xs,{xi_vec,eta_vec});
GRID.y = fnval(OUT.ys,{xi_vec,eta_vec});
% enforce continuity (?)
GRID.x(1:i1,1) = 0.5*( GRID.x(1:i1,1) + GRID.x(end:-1:i2,1) );
GRID.x(end:-1:i2,1) = GRID.x(1:i1,1);

GRID.y(1:i1,1) = 0.5*( GRID.y(1:i1,1) + GRID.y(end:-1:i2,1) );
GRID.y(end:-1:i2,1) = GRID.y(1:i1,1);

% theta = OUT.base_grid.theta;
% nt = numel(theta);
% tfun = @(t) interp1(linspace(0,1,nt),theta,t,"pchip");
tfun = OUT.base_grid.theta_fun;

% need to adjust theta slightly to make the points orthogonal
fz = @(t) airfoil.zs_from_theta(tfun(t));
options = optimset('TolFun',1e-15,'TolX',1e-17);
for i = i1+1:i2-1
    t0 = (i-i1)/(i2-i1);
    t = fzero( @(t) ( real(fz(t)) - GRID.x(i,2) ).*imag(airfoil.unit_normal_cmplx(tfun(t))) ...
             - ( imag(fz(t)) - GRID.y(i,2) ).*real(airfoil.unit_normal_cmplx(tfun(t))), ...
             t0, options);
    GRID.x(i,1) = real(fz(t));
    GRID.y(i,1) = imag(fz(t));
end
end

function P2D_grid_out(GRID,filename)
imax = GRID.imax;
jmax = GRID.jmax;
fid = fopen(filename,'w');
intfmt = '        %4d';
fltfmt = ' %-# 23.16E';
% fprintf(fid,[intfmt,'\n'],1); % this code is just for single block grids
fprintf(fid,[intfmt,intfmt,'\n'],imax,jmax);

count = 0;
for j = 1:jmax
    for i = 1:imax
        fprintf(fid,fltfmt,GRID.x(i,j));
        count = count + 1;
        if count == 2
            fprintf(fid,'\n');
            count = 0;
        end
    end
end

for j = 1:jmax
    for i = 1:imax
        fprintf(fid,fltfmt,GRID.y(i,j));
        count = count + 1;
        if count == 2
            fprintf(fid,'\n');
            count = 0;
        end
    end
end

fclose(fid);

end

function GRID = grid_subset_2D(GRID,small_ind)
GRID.x = GRID.x(small_ind{:});
GRID.y = GRID.y(small_ind{:});
GRID.imax = length(small_ind{1});
GRID.jmax = length(small_ind{2});
end