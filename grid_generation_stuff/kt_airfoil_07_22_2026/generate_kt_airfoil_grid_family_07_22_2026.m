%% KT airfoil test case grid generation (07/22/2026)
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


folder = 'C:\Users\wajordan\Downloads\kt_for_transfer_AR_1_1_TE_10_d-01_ALPHA_5';
out_folder = fullfile(folder,'\grids\');
jobfmt = '_kt%0.4dx%0.4d';

levels = 1:6;
n_levels = length(levels);
airfoil_inputs = struct();
airfoil_inputs.epsilon = 0.1;
airfoil_inputs.kappa   = 0.0;
airfoil_inputs.tau     = 0.0; % deg2rad(10.0);
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = 0.0; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;

n_refine_i   = 16;
n_refine_j   = 16;
n_body_pts = 129;
n_wake_pts = 129;
jmax       = 129;
imax       = n_body_pts + 2*(n_wake_pts-1);
AR_LE      =  1.0;
AR_TE      =  1.0;
delta_LE   = 0.125;
delta_TE   = 0.125;
spline_order = 5;
% these work well for
% jmax == n_body_pts == n_wake_pts == 129
F1 = transfinite_scalar_interpolant_2D( G1 = 1.0000,                        ...
                                        G2 = @(t) smooth_ramps(t,           ...
                                               [ 0.000, 0.2500, 1.0000 ],   ...
                                               [ 1.000, 0.9999, 0.9800 ]),  ...
                                        G3 = 0.9800,                        ...
                                        G4 = @(t) smooth_ramps(t,           ...
                                               [ 0.000, 0.2500, 1.0000 ],   ...
                                               [ 1.000, 0.9999, 0.9800 ]) );
F2 = transfinite_scalar_interpolant_2D( G1 = 0.5,                           ...
                                        G2 = @(t) smooth_ramps(t,           ...
                                               [ 0.00, 0.25, 0.75, 1.00 ],  ...
                                               [ 0.50, 4.00, 4.00, 0.50 ]), ...
                                        G3 = 0.5,                           ...
                                        G4 = @(t) smooth_ramps(t,           ...
                                               [ 0.00, 0.25, 0.75, 1.00 ],  ...
                                               [ 0.50, 4.00, 4.00, 0.50 ]) );
scjmax = @(i,j) F1.eval((i-1)/(imax-1),(j-1)/(jmax-1));
alpham = @(i,j) F2.eval((i-1)/(imax-1),(j-1)/(jmax-1));

i_fine = n_refine_i*(imax-1)+1;
j_fine = n_refine_j*(jmax-1)+1;
[OUT,airfoil,fine_grid] = generate_kt_airfoil_grid_07_22_2026(       ...
                                                     airfoil_inputs, ...
                                                     n_body_pts,     ...
                                                     n_wake_pts,     ...
                                                     jmax,           ...
                                                     scjmax,         ...
                                                     alpham,         ...
                                                     AR_LE,          ...
                                                     AR_TE,          ...
                                                     delta_LE,       ...
                                                     delta_TE,       ...
                                                     spline_order,   ...
                                                     i_fine,        ...
                                                     j_fine );
GRID = fine_grid.GRID;
i1   = fine_grid.i1;
i2   = fine_grid.i2;


dstr = char(datetime('now',"Format",'uuuu-MM-dd''_''HH-mm-ss'));
% save(['OUT_STRUCT_',dstr],"OUT");

skip1 = 2;

hold on
plot( OUT.base_grid.x(1:2^skip1:end,1:2^skip1:end),  ...
      OUT.base_grid.y(1:2^skip1:end,1:2^skip1:end),'k')
plot( OUT.base_grid.x(1:2^skip1:end,1:2^skip1:end).',...
      OUT.base_grid.y(1:2^skip1:end,1:2^skip1:end).','k')
axis equal

skip2 = 8;
hold on;
plot( GRID.x(1:2^skip2:end,1:2^skip2:end), ...
      GRID.y(1:2^skip2:end,1:2^skip2:end),'r');
plot( GRID.x(1:2^skip2:end,1:2^skip2:end).', ...
      GRID.y(1:2^skip2:end,1:2^skip2:end).','r');
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