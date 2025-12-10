%% Generate family of grids for circular obstacle (03/29/23)

clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\svf\svf\';
jobfmt = '_svf%0.4dx%0.4d';

% Nfine = [2049,1025];
Nfine = [1025,513];
N_levels = 8;

GRID = svf_grid(Nfine(1),Nfine(2));
for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset(GRID,{1:s:Nfine(1),1:s:Nfine(2),1});
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\svf'];
    status = mkdir(foldername);
%     P3D_grid_out(GRID2,[filename,'.grd'])
%     svf_bc_out_mms(GRID2.imax,GRID2.jmax,[filename,'.bc'])
    svf_bc_out(GRID2.imax,GRID2.jmax,[filename,'.bc'])
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

function svf_bc_out(N1,N2,filename)
fid = fopen(filename,'w');
intfmt = '        %-5d';
% fltfmt = ' %-# 23.16E';

fprintf(fid,'%d\n',1);
fprintf(fid,'%d\n',1);
fprintf(fid,'%d     %d\n',N1,N2);
fprintf(fid,'A\n');
fprintf(fid,'%d\n',4);
fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1,-200);
fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'],N1,N1, 1,N2, 601);
fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,N1,N2,N2,-200);
fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'], 1, 1, 1,N2,-200);

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