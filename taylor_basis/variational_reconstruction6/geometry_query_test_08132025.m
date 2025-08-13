%% Geometry Query Test (08/13/2025)
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
agglom=true;
load_file=true;
GRID = load_gen_svf_grid_for_testing(parent_dir,agglom,load_file);

dim     = 2;
degree  = 3;
n1      = 1;
n_vars  = 1;

gamma  = 1.4;
inputs = [2.0, 1.0, 0.8611583247416177E+05, 2.0];
ref_inputs = [1.0, 1.0, 347.2206293753677000 ];
test_funs{1} = @(x,y) uvel_svf(x,y,0,gamma,inputs(:),ref_inputs(:));

blk = 1;
idx_low  = [1,1,1];
idx_high = [10,10,1];
SUB_GRID = GRID.subset_grid(blk,idx_low,idx_high);

omega  = 1.3;
n_iter = 200;

C = var_rec_t6_array_helper( SUB_GRID, 1, [1,1,1], SUB_GRID.gblock.Ncells, degree, test_funs, n1, 'vector_dist' );
tic
% C = C.perform_reconstruction_fully_coupled();
C = C.perform_iterative_reconstruction_SOR(omega,n_iter);
t1 = toc

CELLS = var_rec_t5.set_up_cell_var_recs(SUB_GRID,1,[1,1,1],SUB_GRID.gblock.Ncells,degree,test_funs,1,false,'vector_dist');
tic
% CELLS = var_rec_t5.perform_reconstruction_fully_coupled(1,CELLS,[]);
CELLS = var_rec_t5.perform_iterative_reconstruction_SOR(1,CELLS,omega,n_iter);
t2 = toc


% GRID_VARS = grid_vars6(SUB_GRID,blk);
% 
% tol = 1e-13;
% 
% for i = 1:GRID_VARS.ncells
%     face_ids = CELLS(i).nbor_face_id;
%     for j = 1:numel(face_ids)
%         PTS_A = CELLS(i).fquad( CELLS(i).nbor_face_id(j) ).quad_pts;
%         TMP   = GRID_VARS.get_face_quads(GRID,i,true);
%         PTS_B = TMP{j}.quad_pts;
%         diff = abs( PTS_A - PTS_B);
%         if (any(diff > tol))
%             error('mismatching quadrature points!')
%         end
% 
%         PTS_A = CELLS(i).fvec( CELLS(i).nbor_face_id(j) ).v;
%         TMP   = GRID_VARS.get_face_normals(GRID,i,true);
%         PTS_B = TMP{j}.v;
%         diff = abs( PTS_A - PTS_B);
%         if (any(diff > tol))
%             error('mismatching normal vectors!')
%         end
%     end
% end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GRID = load_gen_svf_grid_for_testing(parent_dir,agglom,load_file)

if agglom
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_svf_agglom.mat');
else
    grid_t_file = fullfile(parent_dir,'grid_geometry','grid_derived_type','GRID_svf.mat');
end

if isfile(grid_t_file)&&load_file
    load(grid_t_file,"GRID");
else
    GRID = svf_grid(21,21);
    if agglom
        GRID = grid_type(GRID,agglomerate=true,calc_quads=true,nquad=3,nskip=[2,2,1]);
    else
        GRID = grid_type(GRID,calc_quads=true,nquad=3);
    end
    save(grid_t_file,"GRID");
end

end

function GRID = svf_grid(N_t,N_r)

GRID.dim  = 2;
N_z = 1;
GRID.imax = N_t;
GRID.jmax = N_r;
GRID.kmax = N_z;
t_range = [pi/2, 0];
r_range = [2, 3];
z_range = [0, 0];

t = linspace(t_range(1),t_range(2),N_t);
r = linspace(r_range(1),r_range(2),N_r);
z = linspace(z_range(1),z_range(2),N_z);

[T,R,GRID.z] = ndgrid(t,r,z);

GRID.x = R.*cos(T);
GRID.y = R.*sin(T);

end