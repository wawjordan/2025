%% KT airfoil test case
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\kt_new_grid\no_scale_rot\';
jobfmt = '_kt%0.4dx%0.4d';

inputs = struct();
inputs.l          = 1;
inputs.epsilon    = 0.1;
inputs.kappa      = -0.0;
% inputs.tau        = deg2rad(0.174532925199428);
inputs.tau        = 0;
inputs.r_outer    = 100;
inputs.scale      = true;
inputs.theta_start = 0;
inputs.theta_end = 2*pi;

% inputs.l          = 1;
% inputs.epsilon    = 0.1;
% inputs.kappa      = -0.0;
% inputs.tau        = deg2rad(10);
% inputs.r_outer    = 100;
% inputs.scale      = false;
% inputs.theta_start = pi/2;
% inputs.theta_end = 5*pi/2;

% inputs.l          = 1;
% inputs.epsilon    = 0.0;
% inputs.kappa      = 0.0;
% inputs.tau        = 0.0;
% inputs.r_outer    = 2;
% inputs.scale      = false;
% inputs.theta_start = pi/2;
% inputs.theta_end = 2*pi/2;



% levels = 1:10;
levels = 3:10;
% levels = 1:2;
N_levels = length(levels);
% N_r     = 769;
% N_theta = 1025;
N_r     = 3073;
N_theta = 4097;
% N_theta = 2049;

Nfine = [N_r,N_theta];
[c,~,~] = get_airfoil_chord(inputs.l,inputs.epsilon,inputs.kappa,inputs.tau);
expansion_ratio = (c*inputs.r_outer)^(1/N_r);
inputs.spacing{1} = @(xmin,xmax,N)mygeomspace( xmin, xmax, expansion_ratio, N );
% inputs.spacing{1} = @(xmin,xmax,N)mygeomspace( xmin, xmax, 1.01, N );
inputs.spacing{2} = @linspace;


fprintf('Generating fine grid...\n')
GRID = kt_grid(N_r,N_theta,inputs);

skip = 1;
GRID.x = GRID.x(1:2^(skip-1)+1,:);
GRID.y = GRID.y(1:2^(skip-1)+1,:);
hold on;
plot(GRID.x(1:2^(skip-1):end,1:2^(skip-1):end),GRID.y(1:2^(skip-1):end,1:2^(skip-1):end),'r');
plot(GRID.x(1:2^(skip-1):end,1:2^(skip-1):end).',GRID.y(1:2^(skip-1):end,1:2^(skip-1):end).','r');
axis equal

% inner_wall =  201; % slip wall
% outer_wall = [401   2.0044593143438000e-1  1.0e5  3.4746685537239983e2  0.0  0.0];
% inner_wall = -200;
% outer_wall = -200;

load('C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\kt_new_grid\xy_airfoil_coords_101524.mat')

plot(x_airfoil,y_airfoil,'k.-')


% kt_bc_out(GRID.imax,GRID.jmax,inner_wall,outer_wall,'test.bc')
bc_id_list  = [-1,-1,201,-200];
for j = 1:length(levels)
    s = 2^(levels(j)-1);
    fprintf('writing level %d of %d\n',j,N_levels)
    GRID2 = grid_subset_2D(GRID,{1:s:Nfine(1),1:s:Nfine(2)});
    foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
    filename = [foldername,'\kt'];
    status = mkdir(foldername);
    P2D_grid_out(GRID2,[filename,'.grd'])
%     kt_bc_out(GRID2.imax,GRID2.jmax,inner_wall,outer_wall,[filename,'.bc'])
    N = [GRID2.imax,GRID2.jmax];
    idx_list{1} = [ 1, N(1),  1,  1 ];
    idx_list{2} = [ 1, N(1), N(2), N(2) ];
    idx_list{3} = [ 1,  1,  1, N(2) ];
    idx_list{4} = [N(1), N(1),  1, N(2) ];
    additional_inputs{1} = [ 1, N(1), N(2), N(2),  1 ];
    additional_inputs{2} = [ 1, N(1),  1,  1,  1 ];
    additional_inputs{3} = [];
    additional_inputs{4} = [];
    SENSEI_BC_write(N,idx_list,bc_id_list,additional_inputs,[filename,'.bc'])
end


function GRID = kt_grid(N_r,N_t,inputs)
l       = inputs.l;
epsilon = inputs.epsilon;
kappa   = inputs.kappa;
tau     = inputs.tau;
r_outer = inputs.r_outer;

% radius of circle in zeta plane
a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
n = 2 - tau/pi;
zeta_c = -l*epsilon + l*kappa*1i;
beta = asin(kappa/a0);

[c,zLE,~] = get_airfoil_chord(l,epsilon,kappa,tau);
r     = inputs.spacing{1}(a0,c*r_outer,N_r);
% theta = inputs.spacing{2}(0,2*pi,N_t);
theta = inputs.spacing{2}(inputs.theta_start,inputs.theta_end,N_t);
[R,THETA] = ndgrid(r,theta);
ZETA = R.*exp(1i*(THETA - beta)) + zeta_c;
Z = KT_zeta_to_z(ZETA,n,l);

if inputs.scale
    Z = (Z - zLE)/c;
end

GRID.imax = N_r;
GRID.jmax = N_t;
GRID.x = real(Z);
GRID.y = imag(Z);
end

function [c,zLE,zTE] = get_airfoil_chord(l,epsilon,kappa,tau)
% radius of circle in zeta plane
a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
beta = asin(kappa/a0);
n = 2 - tau/pi;
zeta_c = -l*epsilon + l*kappa*1i;
z_fun = @(s,theta)s*real(KT_zeta_to_z(a0.*exp(1i*(theta - beta)) + zeta_c,n,l));
theta_LE = fminbnd(@(theta)z_fun( 1,theta),0,2*pi);
theta_TE = fminbnd(@(theta)z_fun(-1,theta),0,2*pi);

zeta_LE = a0.*exp(1i*(theta_LE - beta)) + zeta_c;
zeta_TE = a0.*exp(1i*(theta_TE - beta)) + zeta_c;
zLE = KT_zeta_to_z(zeta_LE,n,l);
zTE = KT_zeta_to_z(zeta_TE,n,l);
c = ( abs(zLE) + abs(zTE) );
end

function z = KT_zeta_to_z(zeta,n,l)
zeta_p = (zeta+l).^n;
zeta_m = (zeta-l).^n;
z = n*l*( (zeta_p + zeta_m)./(zeta_p - zeta_m) );
end


function kt_bc_out(N1,N2,inner_wall,outer_wall,filename)

N = [N1,N2];
idx_list{1} = [ 1, N1,  1,  1 ];
idx_list{2} = [ 1, N1, N2, N2 ];
idx_list{3} = [N1, N1,  1, N2 ];
idx_list{4} = [ 1,  1,  1, N2 ];
bc_id_list  = [-1,-1,inner_wall(1),outer_wall(1)];
additional_inputs{1} = [ 1, N1, N2, N2,  1 ];
additional_inputs{2} = [ 1, N1,  1,  1,  1 ];
additional_inputs{3} = inner_wall(2:end);
additional_inputs{4} = outer_wall(2:end);

SENSEI_BC_write(N,idx_list,bc_id_list,additional_inputs,filename)

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

function x = mygeomspace( xmin, xmax, r, N )
% geometric spacing -> each subsequent interval is r times the length of
% previous
% e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
gvec = r.^(0:N-2);
dx = ( xmax - xmin )/sum(gvec);
x = zeros(1,N);
x(1) = xmin;
for i = 2:N
    x(i) = x(i-1) + (r^(i-2)*dx);
end
end

% function kt_bc_out(N1,N2,inner_wall,outer_wall,filename)
% fid = fopen(filename,'w');
% intfmt = '        %-4d';
% fltfmt = ' %-# 23.16E';
% 
% fprintf(fid,'%d\n',1);
% fprintf(fid,'%d\n',1);
% fprintf(fid,'%d     %d\n',N1,N2);
% fprintf(fid,'A\n');
% fprintf(fid,'%d\n',4);
% fprintf(fid,['bc1   ',repmat(intfmt,[1,10]),'\n'], 1,  N1, ...
%                                                   1,   1, ...
%                                                  -1,      ...
%                                                   1,  N1, ...
%                                                   N2, N2, 1 );
% fprintf(fid,['bc2   ',repmat(intfmt,[1,10]),'\n'], 1,  N1, ...
%                                                   N2, N2, ...
%                                                  -1,      ...
%                                                   1,  N1, ...
%                                                   1,   1, 1 );
% fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,   1,  ...
%                                                   1,  N2, ...
%                                                   inner_wall );
% fmt = [repmat(intfmt,[1,5]),repmat(fltfmt,1,length(outer_wall)-1)];
% fprintf(fid,['bc4   ',fmt,'\n'], N1, N1, 1,  N2, outer_wall{:} );
% 
% fclose(fid);
% end