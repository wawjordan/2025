%% Potential Ellipse Testing (10/03/2026)
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

e = potential_ellipse(a=3,b=0.1,alpha=90);
imax = 129;
jmax = 129;
bdist = 10;
% h = e.plot_ellipse();
% axis equal
GRID = e.generate_mapped_grid(bdist,imax,jmax);
% GRID = e.extruded_grid(imax,jmax);
X = GRID.x;
Y = GRID.y;

x = linspace( -e.a/2, e.a/2, imax );
y = linspace( -e.a, 0.0, jmax );
[X,Y] = ndgrid(x,y);
U = e.x_velocity(X,Y);
V = e.y_velocity(X,Y);
P = e.pressure(X,Y);

hold on;
contourf(X,Y,P,'EdgeColor','none')
% quiver(X,Y,U,V)
ns = 101;
ys = repmat(-e.a,1,ns);
xs = linspace(-e.a/2,e.a/2,ns);
verts = stream2(X.',Y.',U.',V.',xs,ys);
streamline(verts,'Color','k');
e.plot_ellipse;
axis equal