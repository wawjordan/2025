%% parabolic stagnation flow testing 10/07/2026
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
in.rLE     = 1.0;
in.rho_inf = 1.0;
in.p_inf   = 100000.0;
in.v_inf   = 68.0;
in.stag_spacing=in.rLE/20;
in.boundary_distance=500;
in.AR = 1;
in.imax=129;
in.jmax=65;

parab = parabolic_stagnation_flow( rLE=in.rLE, ...
                               rho_inf=in.rho_inf, ...
                               v_inf=in.v_inf,     ...
                               p_inf=in.p_inf );

GRID =  parab.extruded_grid( in.imax, ...
                         in.jmax, ...
                         in.stag_spacing, ...
                         in.boundary_distance, ...
                         in.AR );
% GRID2 = parab.parabolic_grid( in.imax, ...
%                           in.jmax, ...
%                           in.stag_spacing, ...
%                           in.boundary_distance, ...
%                           in.AR );

x1 = GRID.x(:,:);
y1 = GRID.y(:,:);

% x2 = GRID2.x(:,:);
% y2 = GRID2.y(:,:);

% hold on;
% plot(x1,y1,'k');
% plot(x1.',y1.','k');
% plot(x2,y2,'r');
% plot(x2.',y2.','r');
% axis equal

u   = parab.x_velocity(GRID.x,GRID.y);
v   = parab.y_velocity(GRID.x,GRID.y);
p   = parab.pressure(GRID.x,GRID.y);
rho = parab.density(GRID.x,GRID.y);


hold on
contourf(GRID.x,GRID.y,p,'EdgeColor','none')
colorbar
plot(x1,y1,'k');
plot(x1.',y1.','k');
axis equal

% quiver(GRID.x,GRID.y,u,v,0.1,"filled",'k')
% ind = 2;
% verts = stream2(GRID.x.',GRID.y.',-u.',-v.',GRID.x(ind,:),GRID.y(ind,:));
% verts = [verts,stream2(GRID.x.',GRID.y.',u.',v.',GRID.x(ind,:),GRID.y(ind,:))];
% verts = [verts,stream2(GRID.x.',GRID.y.',-u.',-v.',GRID.x(:,1),GRID.y(:,1))];
% verts = [verts,stream2(GRID.x.',GRID.y.',-u.',-v.',GRID.x(:,jmax),GRID.y(:,jmax))];
% streamline(verts,'Color','k');