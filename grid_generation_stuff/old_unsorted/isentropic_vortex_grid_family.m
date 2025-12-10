%% Isentropic Vortex Test Case from David Craig's Dissertation
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\isentropic_vortex\iv\';
jobfmt = '_iv%0.4dx%0.4d';

% Nfine = [1025,513];

N_levels = 6;
Nfine = [513,513];
GRID = isvortex_grid(Nfine(1),Nfine(2));
% hold on;
% plot(GRID.x,GRID.y,'k')
% plot(GRID.x.',GRID.y.','k')
% axis equal

for j = 1:N_levels
    s = 2^(j-1);
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\iv'];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
    isvortex_mms_bc_out(GRID2.imax,GRID2.jmax,[filename,'.bc'])
end


function GRID = isvortex_grid(N_i,N_j)

GRID.imax = N_i;
GRID.jmax = N_j;

% (xi,eta)=>(x,y)square linear grid in [0,1]^2
xi1  = linspace(1,0,N_i); % =>theta
eta1 = linspace(0,1,N_j); % =>r
[XI1,ETA1] = ndgrid(xi1,eta1);
XI2  = XI1  + (1/10)*sin(pi*XI1).*sin(pi*ETA1);
ETA2 = ETA1 + (1/10)*exp(1-ETA1).*sin(pi*XI1-(3/4)).*sin(pi*ETA1);


T = (pi/2)*XI2;
R = 2*ETA2 + 1;
GRID.x = R.*cos(T);
GRID.y = R.*sin(T);

end

function isvortex_wall_bc_out(N1,N2,filename)
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

function isvortex_mms_bc_out(N1,N2,filename)
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