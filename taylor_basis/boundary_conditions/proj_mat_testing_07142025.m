%% script for testing your projection matrices (07/14/2025)
clc; clear; close all;
dim = 3;

origin = zeros(dim,1);

normal = rand(dim,1);
vec1   = rand(dim,1);
vec2   = rand(dim,1);

P = get_P_from_vec(normal);
R = eye(dim) - P;
D = R - P;

vec1_n = P*vec1(:);
vec1_t = R*vec1(:);
vec1_r = D*vec1(:);

vec2_n = P*vec2(:);
vec2_t = R*vec2(:);
vec2_r = D*vec2(:);


ax = gca;
hold on;
axis equal
plot_vec(ax,origin,normal,'k')
plot_vec(ax,origin,vec1  ,'k--')
plot_vec(ax,origin,vec1_n,'r--')
plot_vec(ax,origin,vec1_t,'b--')
plot_vec(ax,origin,vec1_r,'g--')

plot_vec(ax,origin,vec2  ,'k:')
plot_vec(ax,origin,vec2_n,'r:')
plot_vec(ax,origin,vec2_t,'b:')
plot_vec(ax,origin,vec2_r,'g:')
axis equal


function P = get_P_from_vec(vec)
P = vec(:)*vec(:).'/dot(vec,vec);
end

function plot_vec(ax,origin,vec,varargin)
pts = [origin(:),origin(:)+vec(:)];
hold on;
if (numel(vec)==2)
    l = plot(ax,pts(1,:),pts(2,:),varargin{:});
    color = get(l,'Color');
    plot(ax,pts(1,1),pts(2,1),'Color',color,'Marker','o')
    plot(ax,pts(1,2),pts(2,2),'Color',color,'Marker','+')
elseif (numel(vec)==3)
    l = plot3(ax,pts(1,:),pts(2,:),pts(3,:),varargin{:});
    color = get(l,'Color');
    plot3(ax,pts(1,1),pts(2,1),pts(3,1),'Color',color,'Marker','o')
    plot3(ax,pts(1,2),pts(2,2),pts(3,2),'Color',color,'Marker','+')
else
    error('only 2D or 3D supported for plotting')
end
end