%% stagflow2 Testing (10/04/2026)
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

s = stagflow2(length=1,alpha=0);
imax = 129;
jmax = 129;

x = linspace( -2*s.length, 2*s.length, imax );
y = linspace( -2*s.length, 2*s.length, jmax );
[X,Y] = ndgrid(x,y);

U = s.x_velocity(X,Y);
V = s.y_velocity(X,Y);
P = s.pressure(X,Y);

hold on;
contourf(X,Y,P,'EdgeColor','none')
ns = 21;
xs = repmat(s.length/10,1,ns);
ys = linspace(-2*s.length,2*s.length,ns);
verts = stream2(X.',Y.',U.',V.',xs,ys);
verts = [verts,stream2(X.',Y.',-U.',-V.',xs,ys)];
verts = [verts,stream2(X.',Y.',U.',V.',-xs,ys)];
verts = [verts,stream2(X.',Y.',-U.',-V.',-xs,ys)];
streamline(verts,'Color','k');
s.plot_plate('k--');
axis equal

plot(0,s.length*cos(s.alpha),'rx')
plot(0,s.length*cos(pi+s.alpha),'g+')

% hold on;
% contourf(Y.',X.',P.','EdgeColor','none')
% ns = 101;
% ys = repmat(-s.length,1,ns);
% xs = linspace(-s.length/2,s.length/2,ns);
% verts = stream2(Y,X,V,U,ys,xs);
% streamline(verts,'Color','k');
% s.plot_plate('k--');
% axis equal