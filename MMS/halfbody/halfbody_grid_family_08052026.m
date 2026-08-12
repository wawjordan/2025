%% halfbody grid generation 08/05/2026
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

% folder = 'C:\Users\Will\Downloads\halfbody_test_AR_1';
folder = 'C:\Users\wajordan\Downloads\halfbody_test_AR_4_skip';
prefix='hb';
out_folder = fullfile(folder,'\grids\');
jobfmt  = ['_',prefix,'%0.4dx%0.4d'];

in         = struct();
in.xstag   = -1;
in.vinf    = 75.0;
in.rhoinf  = 1.0;
in.pinf    = 100000.0;
in.gamma   = 1.4;
in.rho_ref = 1.0;
in.p_ref   = 100000.0;

body = local_halfbody_generator(in);

levels            = 1:6;
n_levels          = length(levels);
boundary_distance = 500;
stag_spacing      = 1e-3;
AR                = 1;
imax              = 4097;
jmax              = 2049;
% imax              = 257;
% jmax              = 129;

skip1 = 4;
skip2 = 1;

GRID = body.extruded_grid2( imax, jmax, stag_spacing, boundary_distance, AR );
Nfine = [GRID.imax,GRID.jmax];
GRID = grid_subset_2D(GRID,{1:skip1:Nfine(1),1:skip2:Nfine(2)});

s1 = levels(end)-1;
s2 = levels(end)-1;
x = GRID.x(1:2^s1:end-30,1:2^s2:end);
y = GRID.y(1:2^s1:end-30,1:2^s2:end);
hold on;
plot( x,   y,  'k');
plot( x.', y.','k');
axis equal
% xlim([-1,2])
% ylim([-1.5,1.5])

% u = body.x_velocity(x,y);
% v = body.y_velocity(x,y);
% vmag = sqrt( u.^2 + v.^2 );
% p = body.pressure(x,y)-body.pinf;
% contourf(x,y,p,100,'LineStyle','none');

% dstr = char(datetime('now',"Format",'uuuu-MM-dd''_''HH-mm-ss'));
% save([prefix,'OUT_STRUCT_',dstr],"OUT");

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

function body = local_halfbody_generator(in)
a_ref   = sqrt(in.gamma*in.p_ref/in.rho_ref);

% nondimensionalize inputs
vinf   = in.vinf/a_ref;
rhoinf = in.rhoinf/in.rho_ref;
pinf   = in.pinf/(in.rho_ref*a_ref^2);
body        = halfbody( in.xstag, vinf );
body.gamma  = in.gamma;
body.rhoinf = rhoinf;
body.pinf   = pinf;
end