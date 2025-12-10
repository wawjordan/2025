%% Generate family of grids for SVF test case (transposed) (06/28/23)

clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\svf\svf_T\';
jobfmt = '_svf%0.4dx%0.4d';

% Nfine = [1025,513];
Nfine = [513,1025];
N_levels = 8;

% Nfine = [9,17];
% N_levels = 1;

GRID = svf_grid_T(Nfine(1),Nfine(2));
for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\svf'];
    status = mkdir(foldername);
%     P3D_grid_out(GRID2,[filename,'.grd'])
    P2D_grid_out(GRID2,[filename,'.grd'])
%     svf_bc_out_mms(GRID2.imax,GRID2.jmax,[filename,'.bc'])
    svf_bc_out_T(GRID2.imax,GRID2.jmax,[filename,'.bc'])
end

function GRID = svf_grid_T(N_r,N_t)

GRID.dim  = 2;
N_z = 1;
GRID.imax = N_r;
GRID.jmax = N_t;
GRID.kmax = N_z;
t_range = [pi/2, 0];
% r_range = [2, 3];
r_range = [3, 2];

t = linspace(t_range(1),t_range(2),N_t);
r = linspace(r_range(1),r_range(2),N_r);

[R,T] = ndgrid(r,t);

GRID.x = R.*cos(T);
GRID.y = R.*sin(T);

end

% function GRID = svf_grid_T(N_r,N_t)
% 
% GRID.dim  = 2;
% N_z = 1;
% GRID.imax = N_r;
% GRID.jmax = N_t;
% GRID.kmax = N_z;
% t_range = [pi/2, 0];
% r_range = [2, 3];
% z_range = [0, 0];
% 
% t = linspace(t_range(1),t_range(2),N_t);
% r = linspace(r_range(1),r_range(2),N_r);
% z = linspace(z_range(1),z_range(2),N_z);
% 
% [T,R,GRID.z] = ndgrid(t,r,z);
% 
% GRID.x = R.*cos(T);
% GRID.y = R.*sin(T);
% 
% end

function svf_bc_out_T(N1,N2,filename)
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

function svf_bc_out_mms(N1,N2,filename)
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