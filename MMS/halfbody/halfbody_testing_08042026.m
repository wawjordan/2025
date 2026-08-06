%% halfbody grid testing 08/04/2026
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

in = struct();
in.xstag = -1;
in.vinf    = 75.0;
in.rhoinf  = 1.0;
in.pinf    = 100000.0;
in.gamma   = 1.4;
in.rho_ref = 1.0;
in.p_ref   = 100000.0;
body = local_halfbody_generator(in);

hold on
tmin=pi/4;

h1 = body.plot_halfbody(tmin);

[x,y]=body.surface_coords(1.97603146838258);
plot(x,y,'r+')
[x,y]=body.surface_coords(4.30715383879700);
plot(x,y,'g+')
% h2 = body.plot_halfbody_normals(linspace(tmin,2*pi-tmin,65));
[Fx,Fy] = body.forces(tmin,2*pi-tmin);
boundary_distance = 50;
stag_spacing = 1e-1;
AR = 1;
imax = 257;
jmax = 129;
n_ref = 0;

% imax = 2049;
% jmax = 1025;
% n_ref = 5;

GRID2 = body.extruded_grid2(imax,jmax,stag_spacing,boundary_distance,AR);
x = GRID2.x(1:2^n_ref:end,1:2^n_ref:end);
y = GRID2.y(1:2^n_ref:end,1:2^n_ref:end);
hold on
u = body.x_velocity(x,y);
v = body.y_velocity(x,y);
vmag = sqrt( u.^2 + v.^2 );
p = body.pressure(x,y)-body.pinf;
contourf(x,y,p,100,'LineStyle','none');
[x,y]=body.surface_coords(1.97603146838258);
plot(x,y,'r+')
[x,y]=body.surface_coords(deg2rad(63));
plot(x,y,'b+')

plot(x,y,'k');
plot(x.',y.','k');
axis equal
xlim([-6,4])
ylim([-5,5])

% stag_spacing = abs(xstag)/20;
% mu = 0.4;
% GRID1 = f.extruded_grid1(imax,jmax,boundary_distance,mu,AR);
% hold on
% plot(GRID1.x(1:2^n_ref:end,1:2^n_ref:end),  GRID1.y(1:2^n_ref:end,1:2^n_ref:end),  'k');
% plot(GRID1.x(1:2^n_ref:end,1:2^n_ref:end).',GRID1.y(1:2^n_ref:end,1:2^n_ref:end).','k');


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