%% rotation matrix from vectors (06/25/2026)
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

dim = 12;

origin = zeros(dim,1);

x   = rand(dim,1);
y   = rand(dim,1);

u   = rand(dim,1);

Rxy = vector2vector(x,y);


% check that it rotated x to y
( y(:)/vmag(y) ) - Rxy*( x(:)/vmag(x) )

% ax = gca;
% hold on;
% axis equal
% plot_vec(ax,origin,x,'r')
% plot_vec(ax,origin,y,'b')
% plot_vec(ax,origin,u,'k--')
% plot_vec(ax,origin,Rxy*u,'k')
% axis equal

function Rxy = vector2vector(x,y)
% https://arxiv.org/pdf/1801.01724
% note: u,v in S^n (unit sphere, so they need to be unit length, v /= -u
x = x(:)/vmag(x);
y = y(:)/vmag(y);
yxt = y*x.';
xyt = x*y.';
xxt = x*x.';
yyt = y*y.';
xyd = dot(x,y);
Rxy = eye(numel(x)) + yxt - xyt + (1/(1+xyd)) * ( xyd*(yxt + xyt) - (xxt + yyt) );
% Kxy = yxt - xyt;
% Rxy = eye(numel(x)) + Kxy + (1/(1+xyd)) * (yxt - xyt)^2;
end

function val = vmag(v)
val = sqrt( sum(v.^2) );
end

function plot_vec(ax,origin,vec,varargin)
if (numel(vec)==2)
    quiver(ax,origin(1),origin(2),vec(1),vec(2),varargin{:});
elseif (numel(vec)==3)
    quiver3(ax,origin(1),origin(2),origin(3),vec(1),vec(2),vec(3),varargin{:});
else
    error('only 2D or 3D supported for plotting')
end
end