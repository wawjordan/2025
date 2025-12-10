%% Script to generate circular cylinder O-grids from Krivodonova & Berger
% (03/31/2023)
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\krivodonova_circ\New\';
jobfmt = 'circ%03dx%03d';

r0 = 0.5;
rinf = 20;
alpha = 1.1648336;

Nfine = [513,129];
N_levels = 6;

% Nfine = [129,33];
% N_levels = 4;

Mach = 0.05;
Pressure = 101325;

Temperature = 300;
gamma = 1.4;
air_mol = 28.97;
R_u     = 8314.4621;
R_gas   = R_u/air_mol;
sound_speed = sqrt(gamma*R_gas*Temperature);
vmag = Mach*sound_speed;
rho  = gamma*Pressure/(sound_speed^2);

GRID = circ_o_grid(Nfine(1),Nfine(2),r0,rinf,alpha);
for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset(GRID,{1:s:Nfine(1),1:s:Nfine(2),1});
%     plot(GRID2.x(:,:,1),GRID2.y(:,:,1),'k')
%     axis equal
%     hold on; plot(GRID2.x(:,:,1).',GRID2.y(:,:,1).','k')
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\circ'];
    status = mkdir(foldername);
    P3D_grid_out(GRID2,[filename,'.grd'])
    circ_o_grid_bc_out(GRID2.imax,GRID2.jmax,Mach,Pressure,Temperature,[filename,'.bc'])
% imax = (GRID2.imax-1)/2 + 1;
% jmax = (GRID2.jmax-1)/2 + 1;
% circ_o_grid_bc_out(imax,jmax,Mach,Pressure,[filename,'.bc'])
end

% function GRID = circ_o_grid(Ni,Nj,r0,rinf,alpha)
% GRID.imax = Ni;
% GRID.jmax = Nj;
% GRID.kmax = 1;
% GRID.dim  = 2;
% N0 = 33;
% 
% t = [0,linspace(2*pi*(1 - 1/(Ni-1)),0,Ni-1)];
% % r = geomspace(r0,rinf,alpha,Nj);
% r = @(N) interp1( linspace(r0,rinf,N0), geomspace(r0,rinf,alpha,N0), ...
%                   linspace(r0,rinf,N), "pchip" );
% 
% [T,R] = ndgrid(t,r(Nj));
% 
% 
% GRID.x = R.*cos(T);
% GRID.y = R.*sin(T);
% GRID.z = zeros(Ni,Nj);
% end

function GRID = circ_o_grid(Ni,Nj,r0,rinf,alpha)
GRID.imax = Ni;
GRID.jmax = Nj;
GRID.kmax = 1;
GRID.dim  = 2;

t = [0,linspace(2*pi*(1 - 1/(Ni-1)),0,Ni-1)];
r = geomspace(r0,rinf,alpha,Nj);

[T,R] = ndgrid(t,r);


GRID.x = R.*cos(T);
GRID.y = R.*sin(T);
GRID.z = zeros(Ni,Nj);
end

function circ_o_grid_bc_out(N1,N2,Mach,Pressure,Temperature,filename)
fid = fopen(filename,'w');
intfmt = '        %-5d';
fltfmt = ' %-# 23.16E';

fprintf(fid,'%d\n',1);
fprintf(fid,'%d\n',1);
fprintf(fid,'%d     %d\n',N1,N2);
fprintf(fid,'A\n');
fprintf(fid,'%d\n',4);
fprintf(fid,['connected1',repmat(intfmt,[1,10]),'\n'], 1, 1, 1,N2,-1,N1,N1, 1,N2, 1);
fprintf(fid,['connected2',repmat(intfmt,[1,10]),'\n'],N1,N1, 1,N2,-1, 1, 1, 1,N2, 1);
fprintf(fid,['slipwall  ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1, 201);
fprintf(fid,['farfield  ',[repmat(intfmt,[1,5]),repmat(fltfmt,[1,6])],'\n'], 1,N1,N2,N2, 401, Mach, Pressure, Temperature, 0, 0);
fclose(fid);
end

% r1 = geomspace(0.5,20,alpha,33);
% r2 = circ_rad(0.5,33,alpha);
% function r = circ_rad(r0,N1,alpha)
% 
% r = zeros(1,N1);
% r(1) = r0;
% for j = 2:N1
%     tmp = 0;
%     for k = 0:j-2
%         tmp = tmp + alpha^k;
%     end
%     r(j) = r0*(1 + (2*pi/128)*tmp);
% end
% 
% end