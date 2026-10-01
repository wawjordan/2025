%% stagnation flow testing 09/30/2026
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
in.x0      = 0;
in.y0      = 0;
in.x_max   = -1.0;
in.y_max   =  0.5;
in.a_ref   = 300;
in.rho_ref = 1.0;
in.rho_inf = 1.0;
in.gamma   = 1.4;
in.p_stag  = 100000.0;
in.p_ref   = 100000.0;
in.p_min   =  80000.0;
in.M_max   = 0.5;
stag = stagflow( x0=in.x0, ...
                 y0=in.y0, ...
                 lx=in.x_max,...
                 ly=in.y_max,...
                 aref=in.a_ref,...
                 rhoref=in.rho_ref,...
                 rhoinf=in.rho_inf,...
                 gamma=in.gamma,...
                 pstag=in.p_stag,...
                 pmin=in.p_min,...
                 mach_max=in.M_max);

GRID = stag.make_grid(65,33);

u   = stag.x_velocity(GRID.x,GRID.y);
v   = stag.y_velocity(GRID.x,GRID.y);
p   = stag.pressure(GRID.x,GRID.y);
rho = stag.density(GRID.x,GRID.y);
m   = stag.mach(GRID.x,GRID.y);


ind = 2;
hold on
contourf(GRID.x,GRID.y,m,'EdgeColor','none')
% quiver(GRID.x,GRID.y,u,v,1,"filled",'k')
verts = stream2(GRID.x.',GRID.y.',-u.',-v.',GRID.x(ind,:),GRID.y(ind,:));
verts = [verts,stream2(GRID.x.',GRID.y.',u.',v.',GRID.x(ind,:),GRID.y(ind,:))];
verts = [verts,stream2(GRID.x.',GRID.y.',-u.',-v.',GRID.x(:,1),GRID.y(:,1))];
verts = [verts,stream2(GRID.x.',GRID.y.',-u.',-v.',GRID.x(:,end),GRID.y(:,end))];
streamline(verts,'Color','k');
axis equal
colorbar