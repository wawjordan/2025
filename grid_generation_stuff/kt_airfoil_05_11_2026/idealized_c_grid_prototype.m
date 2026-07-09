%% prototype 'idealized' C-grid (06/16/2026)
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

t = linspace(0,1,1001);

% [val,dval,ddval] = combined(t, t0,  t1,  s,  s2,  N,    r )
% [val,dval,ddval] = combined(  t,  0, 0.2, 1, 1.2, 33, 1.05 );

% plot(t,dval)
% plot(t,ddval)

r0 = 1;
r1 = 500+r0;
l0 = 1;

[x1,y1] = body_pts(t,r0,l0);
[x2,y2] = body_pts(t,r1,l0);
x = [x1;x2];
y = [y1;y2];
hold on
plot(x,y,'k')
plot(x.',y.','k')
hold off

% function [x,y] = idealized_c_grid(wall_spacing,r0,l0,r1,l1,N_body,N_wake,jmax)
% end

% function [x,y] = body_pts(t,r,l)
% arc_length = pi*r + 2*l;
% tc1 = l/arc_length;
% tc2 = 1 - tc1;
% x = zeros(size(t));
% y = zeros(size(t));
% mask1 = t<tc1;
% mask2 = (t>=tc1)&(t<=tc2);
% mask3 = t>tc2;
% 
% x(mask1) = l - l*t(mask1)/tc1;
% y(mask1) = -r;
% 
% x(mask2) = r*cos(pi*(0.5+(tc2-t(mask2))/(tc2-tc1)) );
% y(mask2) = r*sin(pi*(0.5+(tc2-t(mask2))/(tc2-tc1)) );
% 
% x(mask3) = l*(t(mask3)-tc2)/(1-tc2);
% y(mask3) =  r;
% end

function [x,y] = body_pts(t,r,l)
arc_length = pi*r + 2*l;
tc1 = l/arc_length;
tc2 = 1 - tc1;
x = zeros(size(t));
y = zeros(size(t));
mask1 = t<tc1;
mask2 = (t>=tc1)&(t<=tc2);
mask3 = t>tc2;

x(mask1) = l - l*t(mask1)/tc1;
y(mask1) = -r;

x(mask2) = r*cos(pi*(0.5+(tc2-t(mask2))/(tc2-tc1)) );
y(mask2) = r*sin(pi*(0.5+(tc2-t(mask2))/(tc2-tc1)) );

x(mask3) = l*(t(mask3)-tc2)/(1-tc2);
y(mask3) =  r;
end



% [x,] = my_geomspace(N,xmin,varargin)
% 
% 
% function [x,r,dx1] = geomspace_4( N, xmin, xmax, dx0 )
% % geometric spacing -> each subsequent interval is r times the length of
% % previous
% % e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
% 
% S = xmax - xmin;
% 
% if S/(N-2) > dx0
%     r0 = 1.01;
% else
%     r0 = 0.99;
% end
% 
% fun = @(r) S/dx0 - sum(r.^(0:N-2));
% options = optimset('FunValCheck','on');
% r = fzero(fun,r0,options);
% 
% x = zeros(1,N);
% x(1) = xmin;
% for i = 2:N-1
%     x(i) = x(i-1) + (r^(i-2)*dx0);
% end
% x(N) = xmax;
% dx1 = x(N)-x(N-1);
% end


% function [x,y] = wake_pts_inner(t,r0,r1,l0,l1)
% arc_length = pi*r0 + 2*l0;
% tc1 = l/arc_length;
% dxdt = (l0)

function [val,dval,ddval] = combined(t,t0,t1,s,s2,N,r)
val   = zeros(size(t));
dval  = zeros(size(t));
ddval = zeros(size(t));

% find x0,dx0
% log(rN)*dx0/(r-1)=s
% dx0 = s*(r-1)/log(r^(N-1));
dx0 = s*(r-1)/log(r^(N-1));
x0  = linear(t2,t0,s);
mask1 = (t<=t0);
mask2 = (t>t0)&(t<=t1);
mask3 = (t>t1)&(t<=t2);
mask4 = (t>t2);
val(mask1) = 0;


[f1,df1,~]    = linear(t1,t0,s);
[f2,df2,ddf2] = stretch(t2,t2,N,r,x0,dx0);
[val(mask2),dval(mask2),ddval(mask2)] = linear(t(mask2),t0,s);
[val(mask3),dval(mask3),ddval(mask3)] = blend(t(mask3),t1,t2,f1,f2,df1,df2,ddf2);
[val(mask4),dval(mask4),ddval(mask4)] = stretch(t(mask4),t2,N,r,x0,dx0);
end

function [val,dval,ddval] = linear(t,t0,s)
val   = s*(t-t0);
dval  = s*ones(size(t));
ddval =  zeros(size(t));
end

function [val,dval,ddval] = blend(t,t1,t2,f1,f2,df1,df2,ddf2)
% assumes t in [t1,t2]
    [val,dval,ddval] = qh_interp(t,t1,t2,f1,f2,df1,df2,0,ddf2);
end

% log(rN)*dx0/(r-1)=s

function [val,dval,ddval] = stretch(t,t2,N,r,x0,dx0)
% assumes t > t2
rN    = r^(N-1);
dxr   = dx0/(r-1);
val   = x0 + dxr*(rN.^(t-t2) - 1);
dval  = dxr*log(rN)*rN.^(t-t2);
ddval = dxr*log(rN)*log(rN)*rN.^(t-t2);
end