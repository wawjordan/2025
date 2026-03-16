%% Getting number of stencil cells
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 40;


clc;
GRID = grid_type(cart_grid(n,n,n));

fn = @(n) arrayfun(@(n) get_n_stencil(GRID,n),n);

fn(3)
fn(10)
fn(26)
fn(64)
fn(130)
fn(232)



function n_stencil_out = get_n_stencil(GRID,n_stencil)
block_info_list(1)           = block_info_t();
block_info_list(1).block_id  = 1;
block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;
print_iter = false;
blk = 1;
idx = floor(GRID.gblock(1).Ncells/2);
[ ~, n_stencil_out, ~ ] = cell_t.build_stencil( blk, idx, n_stencil, GRID, block_info_list, true, print_iter );
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