%% Generate family of Curvilinear grids (09/20/2023)
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\curvilinear_mms\cts\';
jobfmt = '_cts%0.4dx%0.4d';

% Nfine = [513,513];
Nfine = [65,65];
N_levels = 7;

GRID = curv_grid(Nfine(1),Nfine(2));
for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\cts'];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
    bc_out_mms([GRID2.imax,GRID2.jmax],[filename,'.bc'])
end

function bc_out_mms(N,filename)
fid = fopen(filename,'w');
intfmt = '        %-5d';
% fltfmt = ' %-# 23.16E';
N1 = N(1);
N2 = N(2);
if numel(N) == 2
    fprintf(fid,'%d\n',1);
    fprintf(fid,'%d\n',1);
    fprintf(fid,'%d     %d\n',N1,N2);
    fprintf(fid,'A\n');
    fprintf(fid,'%d\n',4);
    fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'],  1, N1,  1,  1, -200);
    fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'], N1, N1,  1, N2, -200);
    fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'],  1, N1, N2, N2, -200);
    fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'],  1,  1,  1, N2, -200);
elseif numel(N) == 3
    N3 = N(3);
    fprintf(fid,'%d\n',1);
    fprintf(fid,'%d\n',1);
    fprintf(fid,'%d     %d     %d\n',N1,N2,N3);
    fprintf(fid,'A\n');
    fprintf(fid,'%d\n',6);
    fprintf(fid,['bc1   ',repmat(intfmt,[1,7]),'\n'],  1,  1,  1, N2,  1, N3, -200);
    fprintf(fid,['bc2   ',repmat(intfmt,[1,7]),'\n'], N1, N1,  1, N2,  1, N3, -200);
    fprintf(fid,['bc3   ',repmat(intfmt,[1,7]),'\n'],  1, N1,  1,  1,  1, N3, -200);
    fprintf(fid,['bc4   ',repmat(intfmt,[1,7]),'\n'],  1, N1, N2, N2,  1, N3, -200);
    fprintf(fid,['bc5   ',repmat(intfmt,[1,7]),'\n'],  1, N1,  1, N2,  1,  1, -200);
    fprintf(fid,['bc6   ',repmat(intfmt,[1,7]),'\n'],  1, N1,  1, N2, N3, N3, -200);
end
fclose(fid);
end

function GRID = curv_grid(N_x,N_y,N_z)

GRID.imax = N_x;
GRID.jmax = N_y;
x_range = [0, 1];
y_range = [0, 1];
GRID.dim  = 2;
   
xi  = linspace(x_range(1),x_range(2),N_x);
eta = linspace(y_range(1),y_range(2),N_y);
zeta = 0;

if nargin == 3
    GRID.kmax = N_z;
    GRID.dim = 3;
    z_range = [0, 1];
    zeta = linspace(z_range(1),z_range(2),N_z);   
end

[XI,ETA,ZETA] = ndgrid(xi,eta,zeta);

GRID.x = xFunc3D(XI,ETA,ZETA);
GRID.y = yFunc3D(XI,ETA,ZETA);

if nargin == 3
    GRID.z = zFunc3D(XI,ETA,ZETA);
end

end

function x = xFunc3D(xi,eta,zeta)
x =  0.10*sin( (2-xi)*pi .*  eta ) ...
  +  0.05*sin(   2   *pi .* zeta ) ...
  +  0.10*eta + 0.05*zeta + 1.00*xi;
end

function y = yFunc3D(xi,eta,zeta)
y = -0.10*sin( (2+eta)*pi .* xi   ) ...
  +  0.05*sin(   2    *pi .* zeta ) ...
  -  0.10*xi + 0.05*zeta + 1.00*eta;
end

function z = zFunc3D(xi,eta,zeta)
z = -0.03*sin(  2    *pi .* xi .* eta ) ...
  -  0.10*xi + 0.05*eta + 1.00*zeta;
end