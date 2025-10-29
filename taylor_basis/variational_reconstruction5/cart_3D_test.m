%% 3D Test (09/15/2025)
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
load_file=false;
GRID = load_gen_grid_for_testing(parent_dir,load_file);

dim     = 2;
degree  = 3;
n1      = 1;
n_vars  = 1;
% test_funs{1}  = @(x,y,z) 999.0 * x - 888.0 * y + 777.0 * z - 666.0;
test_funs{1}  = @(x,y) 999.0 * x - 888.0 * y - 666.0;
blk = 1;
idx_low  = [1,1,1];
idx_high = [4,4,1];
SUB_GRID = GRID.subset_grid(blk,idx_low,idx_high);

CELLS = var_rec_t5.set_up_cell_var_recs(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,1,false,'vector_dist');

omega  = 1.3;
n_iter = 10;
CELLS = var_rec_t5.perform_iterative_reconstruction_SOR(1,CELLS,omega,n_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GRID = load_gen_grid_for_testing(parent_dir,load_file)
grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_cart_3D.mat');

if isfile(grid_t_file)&&load_file
    load(grid_t_file,"GRID");
else
    GRID = cart_grid(5,5,2);
    GRID = grid_type(GRID,agglom=true,calc_quads=true,nquad=3,nskip=[1,1,1]);
    save(grid_t_file,"GRID");
end

end

function GRID = cart_grid(N_x,N_y,N_z)
GRID.dim  = 3;
GRID.imax = N_x;
GRID.jmax = N_y;
GRID.kmax = N_z;
[GRID.x,GRID.y,GRID.z] = ndgrid( linspace(0,1,N_x), ...
                                 linspace(0,1,N_y), ...
                                 linspace(0,1,N_z) );
end