%% Testing Function (03/16/2026)
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

block_num = 1;
n_dim     = 2;
degree    = 4;
n_vars    = 1;
n_quad    = degree;

% load(fullfile(parent_dir,'taylor_basis','new_reconstruct_prototype_03_2026','GRID_svf_agglom.mat'));
n_dim = 1;
GRID.x = my_geomspace(33,0,1,r=1.0);
GRID = grid_type(GRID,calc_quads=true,nquad=n_quad);


% gamma  = 1.4;
% inputs = [2.0, 1.0, 0.8611583247416177E+05, 2.0];
% ref_inputs = [1.0, 1.0, 347.2206293753677000 ];
% test_funs{1}  = @(x,y) dens_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% test_funs{2} = @(x,y) uvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% test_funs{3} = @(x,y) vvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% test_funs{4} = @(x,y) pres_svf(x,y,0,gamma,inputs(:),ref_inputs(:));
% ext_fun = @(x) test_funs{2}(x(1),x(2));

for i = 1:n_vars
    [test_funs{i},~] = generate_random_poly_fun(n_dim,degree);
end
if (n_dim == 1)
    ext_fun = @(x) test_funs{1}(x(1));
elseif (n_dim == 2)
    ext_fun = @(x) test_funs{1}(x(1),x(2));
end

ext_fun = @(x) 1.0 + x - x.^2 - x^3 - x^4;

REC = rec_t( GRID, n_dim, degree, n_vars );

REC = REC.solve_rec( GRID, ext_fun );
% RB = rec_block_t( GRID, block_num, n_dim, degree, n_vars );
% RB = RB.init_cells(RB.p.n_terms,RB.n_vars,1:RB.n_vars);
