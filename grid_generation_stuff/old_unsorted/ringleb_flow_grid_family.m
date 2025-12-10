%% Ringleb Flow test case
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\ringleb_flow\ring\';
jobfmt = '_ring%0.4dx%0.4d';

inputs = struct();
inputs.psi_min = 0.7;
inputs.psi_max = 1.5;
inputs.V_min   = 0.4;

% Nfine = [1025,513];

N_levels = 6;
Nfine = [257,1025];
% Nfine = [17,65];
GRID = ringleb_grid(Nfine(1),Nfine(2),inputs);
tol = 1e-6;
max_iter = 1000;
[GRID.x,GRID.y,R] = winslow_eqns_fd(GRID.x,GRID.y,max_iter,tol);
% hold on;
% plot(GRID.x,GRID.y,'k')
% plot(GRID.x.',GRID.y.','k')
% axis equal

for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\ring'];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
    ringleb_wall_bc_out(GRID2.imax,GRID2.jmax,[filename,'.bc'])
end


function GRID = ringleb_grid(N_i,N_j,inputs)

GRID.imax = N_i;
GRID.jmax = N_j;
[GRID.x,GRID.y,~,~] = ringleb_transfinite_spline_grid_2( N_i, N_j, ...
                                                         inputs.psi_min, ...
                                                         inputs.psi_max, ...
                                                         inputs.V_min );
end

function ringleb_wall_bc_out(N1,N2,filename)
fid = fopen(filename,'w');
intfmt = '        %-5d';
% fltfmt = ' %-# 23.16E';

fprintf(fid,'%d\n',1);
fprintf(fid,'%d\n',1);
fprintf(fid,'%d     %d\n',N1,N2);
fprintf(fid,'A\n');
fprintf(fid,'%d\n',4);
fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1,-200);
fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'],N1,N1, 1,N2, 201);
fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,N1,N2,N2,-200);
fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'], 1, 1, 1,N2, 201);

fclose(fid);
end

function ringleb_mms_bc_out(N1,N2,filename)
fid = fopen(filename,'w');
intfmt = '        %-5d';
% fltfmt = ' %-# 23.16E';

fprintf(fid,'%d\n',1);
fprintf(fid,'%d\n',1);
fprintf(fid,'%d     %d\n',N1,N2);
fprintf(fid,'A\n');
fprintf(fid,'%d\n',4);
fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1,-200);
fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'],N1,N1, 1,N2,-200);
fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,N1,N2,N2,-200);
fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'], 1, 1, 1,N2,-200);

fclose(fid);
end