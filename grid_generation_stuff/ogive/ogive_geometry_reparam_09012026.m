%% ogive geometry reparam work (09/01/2026)
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

rho_ogv = 2.80028;
L0      = 0.5;
Rn      = 12.5/1000;
Rb      = 45/1000;
% y(x) for the ogive part of the geometry
y_ogive    = @(x)   sqrt(rho_ogv^2 - (L0-x).^2) + Rb - rho_ogv;

% dy/dx for the ogive part of the geometry
% dydx_ogive = @(x)   (L0-x)./sqrt(rho_ogv^2 - (L0 - x).^2);

% y(x) (positive branch) for the tip
y_circ     = @(x,x0)   sqrt( Rn^2 - (x-x0).^2 );

% dy/dx for the tip
% dydx_circ  = @(x,x0) -(x-x0)./sqrt(Rn^2- (x-x0).^2);



% center of circle such that it is tangent to the ogive
x0 = L0 - sqrt( (rho_ogv-Rn)^2 - (Rb-rho_ogv)^2 );

xa = x0-Rn;

% point of tangency between circular tip and ogive
options = optimset('TolFun',1e-15,'TolX',1e-17);
xt = fminbnd( @(x) ( (y_circ(x,x0) -  y_ogive(x)).^2 ),x0-Rn,x0+Rn,options);

% either or
yt = y_circ(xt,x0);
% yt = y_ogive(xt);

yfun = @(x) (x<=xt).*y_circ(x,x0) + (x>xt).*y_ogive(x);


hold on
% fplot( @(x) y_circ(x,x0),[x0-Rn,xt],'r')
% fplot(y_ogive,[xt,0.5],'r')
% plot(xt,yt,'gx')
x = linspace(xa,L0,8001);
y = yfun(x);
plot(x,y,'k.-')
axis equal

load('xy_ogive.mat')

plot(xy_ogive(:,1),xy_ogive(:,2),'g');