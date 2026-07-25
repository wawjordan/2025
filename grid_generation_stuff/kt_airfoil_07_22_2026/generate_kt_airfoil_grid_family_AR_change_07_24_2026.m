%% KT airfoil test case grid coarsening for different AR (07/24/2026)
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

folder = 'C:\Users\wajordan\Downloads\kt_AR_1_1_d-1-8_ALPHA_05';
fine_grid_file = 'C:\Users\wajordan\Downloads\kt_AR_1_1_d-1-8_ALPHA_05\grids\_kt6145x2049\kt.grd';

% folder =         'C:\Users\wajordan\Downloads\kt_AR_1_1_d-1-8_ALPHA_00';
% fine_grid_file = 'C:\Users\wajordan\Downloads\kt_AR_1_1_d-1-8_ALPHA_00\grids2\kt_mod.x';

fine_grid = read_grd_file_to_struct(fine_grid_file);

skip1 = 2;
skip2 = 1;

out_folder = fullfile(folder,'\grids_AR_1_2\');
jobfmt = '_kt%0.4dx%0.4d';

levels = 1:6;
n_levels = length(levels);
n_refine_i   = 16;
n_refine_j   = 16;
i1 = 2049;
i2 = 4097;
i1 = (i1-1)/skip1 + 1;
i2 = (i2-1)/skip1 + 1;
GRID = fine_grid.gblock(1);
Nfine = [GRID.imax,GRID.jmax];
GRID = grid_subset_2D(GRID,{1:skip1:Nfine(1),1:skip2:Nfine(2)});

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