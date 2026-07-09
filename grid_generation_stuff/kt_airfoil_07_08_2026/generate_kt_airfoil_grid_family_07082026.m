%% KT airfoil test case (07/08/2026)
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

folder = 'C:\Users\Will\Downloads\kt_for_transfer_07082026_AR_LE_TE_10';
out_folder = fullfile(folder,'\grids\');
jobfmt = '_kt%0.4dx%0.4d';

levels = [1:6];
n_levels = length(levels);
airfoil_inputs = struct();
airfoil_inputs.epsilon = 0.1;
airfoil_inputs.kappa   = 0.0;
airfoil_inputs.tau     = 0.0;
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = 0; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;

n_refine   = 2^4;
n_body_pts = 129;
n_wake_pts = 65;
jmax_in    = 129;
AR_LE      =  10.0;
AR_TE      =  10.0;
delta_LE   = 0.1;
delta_TE   = 0.1;
spline_order = 5;

[OUT,airfoil] = generate_kt_airfoil_grid(airfoil_inputs,n_body_pts,n_wake_pts,jmax_in,AR_LE,AR_TE,delta_LE,delta_TE,spline_order);
hold on
plot(OUT.base_grid.x,OUT.base_grid.y,'k')
plot(OUT.base_grid.x.',OUT.base_grid.y.','k')
axis equal

imax = n_refine*(OUT.base_grid.imax-1)+1;
jmax = n_refine*(OUT.base_grid.jmax-1)+1;
% imax = 311;
% jmax = 111;
[GRID,i1,i2] = get_fine_grid(OUT,airfoil,imax,jmax);
skip1 = levels(end) - 1;
skip2 = max(0,skip1-2);
hold on;
plot(GRID.x(1:2^skip2:end,1:2^skip1:end),GRID.y(1:2^skip2:end,1:2^skip1:end),'r');
plot(GRID.x(1:2^skip1:end,1:2^skip2:end).',GRID.y(1:2^skip1:end,1:2^skip2:end).','r');
axis equal
xlim([-1,2])
ylim([-1.5,1.5])


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
theta = OUT.base_grid.theta;
nt = numel(theta);
tfun = @(t) interp1(linspace(0,1,nt),theta,t,"pchip");
fz = @(t) ( airfoil.airfoil_coords(tfun(t)) - airfoil.airfoil_coords(airfoil.thetaLE) )/airfoil.chord;
for i = i1:i2
    t = (i-i1)/(i2-i1);
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