%% MS1 (inviscid wall MS from Navah & Nadarajah) (08/10/2026)
% https://doi.org/10.1016/j.compfluid.2020.104504
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

folder = 'C:\Users\Will\Downloads\ms1_grids';
% folder = 'C:\Users\wajordan\Downloads\ms1_grids';
prefix='ms1';
out_folder = fullfile(folder,'\grids\');
jobfmt  = ['_',prefix,'%0.4dx%0.4d'];

levels   = 1:6;
n_levels = length(levels);
imax     = 1025;
jmax     = 257;
r_xi     = 0.9993; % 0.993 in paper - but that might be a typo because
r_eta    = 1.0;
a        = 0.05;
N_skip   = 4;      % this will generate grids 4x finer than in the paper, 
                   % but they should comprise the same grid family
GRID     = ms1_grid(imax,jmax,N_skip,r_xi,r_eta,a);
GRID2    = ms1_grid(imax,jmax,1,r_xi,r_eta,a);

skip1 = levels(end)-1;
skip2 = levels(end)-1;
x = GRID.x(1:2^skip1:end,1:2^skip2:end);
y = GRID.y(1:2^skip1:end,1:2^skip2:end);
hold on;
plot( x,   y,  'k');
plot( x.', y.','k');
axis equal

x = GRID2.x(1:2^skip1:end,1:2^skip2:end);
y = GRID2.y(1:2^skip1:end,1:2^skip2:end);
plot( x,   y,  'r');
plot( x.', y.','r');

bc_id_list  = [201,-200,-200,-200];

Nfine = [GRID.imax,GRID.jmax];
for j = 1:n_levels
    s = 2^(levels(j)-1);
    fprintf('writing level %d of %d\n',j,n_levels)
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [out_folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\',prefix];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
    Ni(1) = 1;
    Ni(2) = (GRID.imax - 1)/s + 1;

    Nj(1) = 1;
    Nj(2) = (GRID.jmax - 1)/s + 1;

    N = [GRID2.imax,GRID2.jmax];

    idx_list{1} = [ Ni(1), Ni(2),  Nj(1), Nj(1) ];
    idx_list{2} = [ Ni(1), Ni(2),  Nj(2), Nj(2) ];
    idx_list{3} = [ Ni(1), Ni(1),  Nj(1), Nj(2) ];
    idx_list{4} = [ Ni(2), Ni(2),  Nj(1), Nj(2) ];
    additional_inputs{1} = [];
    additional_inputs{2} = [];
    additional_inputs{3} = [];
    additional_inputs{4} = [];
    SENSEI_BC_write(N,idx_list,bc_id_list,additional_inputs,[filename,'.bc'])
end


function GRID = ms1_grid(imax,jmax,N_refine,r_xi,r_eta,a)
% undeformed initial domain
xi_min = 0.5;
xi_max = 1.0;
eta_min = 0.0;
eta_max = 0.5;

% [xi0,~]  = my_geomspace(imax,xi_min,xmax=xi_max,r=r_xi);
% [eta0,~] = my_geomspace(jmax,eta_min,xmax=eta_max,r=r_eta);

[xi0,~]  = my_geomspace_w_refine(imax,N_refine,xi_min,xmax=xi_max,r=r_xi);
[eta0,~] = my_geomspace_w_refine(jmax,N_refine,eta_min,xmax=eta_max,r=r_eta);

[XI,ETA]  =  ndgrid( xi0, eta0 );

GRID = struct();
GRID.imax = imax;
GRID.jmax = jmax;

GRID.x = x_map(XI,ETA);
GRID.y = y_map(XI,ETA,a);

    function x = x_map(xi,eta)
        x = (4/3)*(xi.^2 -1/4) + 1 + 0*eta;
    end
    function y = y_map(xi,eta,a)
        x = x_map(xi,eta);
        y = eta + a*sin(2*pi*x);
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