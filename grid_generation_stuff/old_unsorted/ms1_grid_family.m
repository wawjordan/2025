%% Isentropic Vortex Test Case from David Craig's Dissertation
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\ms1\ms1\';
jobfmt = '_ms1%0.4dx%0.4d';

Nfine = [1024,256] + 1;

N_levels = 6;
for j = 1:N_levels
    h = 2^(j-1);
    GRID = ms1_grid(h,Nfine);
    foldername = [folder,sprintf(jobfmt,GRID.imax,GRID.jmax)];
    filename = [foldername,'\ms1'];
    status = mkdir(foldername);
    P2D_grid_out(GRID,[filename,'.grd'])
    ms1_wall_bc_out(GRID.imax,GRID.jmax,[filename,'.bc'])
end


function GRID = ms1_grid(h,Nfine)
% undeformed initial domain
xi_range  = [0.5,1.0];
eta_range = [0.0,0.5];

% expansion ratios
r_xi  = 0.9993; % 0.993 in paper - but that must be a typo
r_eta = 1.0;

xi0  = geomspace(  xi_range(1),  xi_range(2),  r_xi,  Nfine(1) );
eta0 = geomspace( eta_range(1), eta_range(2), r_eta,  Nfine(2) );

% coarsened grid
[XI,ETA]  =  ndgrid( xi0(1:h:Nfine(1)), eta0(1:h:Nfine(2)) );

GRID = struct();
GRID.imax = size(XI,1);
GRID.jmax = size(XI,2);

GRID.x = x_map(XI,ETA);
GRID.y = y_map(XI,ETA);
end

function val = a();     val = 0.05; end

function x = x_map(xi,eta)
x = (4/3)*(xi.^2 -1/4) + 1 + 0*eta;
end

function y = y_map(xi,eta)
x = x_map(xi,eta);
y = eta + a()*sin(2*pi*x);
end

function ms1_wall_bc_out(N1,N2,filename)
fid = fopen(filename,'w');
intfmt = '        %-5d';
% fltfmt = ' %-# 23.16E';

fprintf(fid,'%d\n',1);
fprintf(fid,'%d\n',1);
fprintf(fid,'%d     %d\n',N1,N2);
fprintf(fid,'A\n');
fprintf(fid,'%d\n',4);
fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1, 201);
fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'],N1,N1, 1,N2,-200);
fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,N1,N2,N2,-200);
fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'], 1, 1, 1,N2,-200);

fclose(fid);
end