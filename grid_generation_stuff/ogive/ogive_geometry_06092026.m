%% ogive geometry work (06/09/2026)
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
syms rho_ogv positive real
syms L0      positive real
syms Rb      positive real
syms Rn      positive real
syms x       positive real
syms x0      positive real

y_ogv  = sqrt(rho_ogv^2 - (L0-x).^2) + Rb - rho_ogv;
y_circ = sqrt( Rn^2 - (x-x0).^2 );

eqn1 = simplify(y_ogv==y_circ);
eqn1 = simplify(subs(eqn1,{L0,Rb,Rn,rho_ogv},{0.5,45/1000,12.5/1000,2.8}));
s1 = solve(eqn1,x0,'ReturnConditions',true,'Real',true);
for i = 1:numel(s1.x0)
    s1.x0(i) = simplify(s1.x0(i));
end
for i = 1:numel(s1.conditions)
    s1.conditions(i) = simplify(s1.conditions(i));
end


s0 = x - sqrt(Rn^2-(sqrt(rho_ogv^2-(L0-x).^2) + Rb - rho_ogv));
s0 = simplify(subs(s0,{L0,Rb,Rn,rho_ogv},{0.5,45/1000,12.5/1000,2.8}));
f0 =  matlabFunction(s0);

f1 =  matlabFunction(s1.x0(1));
f2 =  matlabFunction(s1.x0(2));

% hold on
% fplot(f1,'r')
% fplot(f2,'b')

dydx_circ = simplify(diff(y_circ,x));
dydx_ogv  = simplify(diff(y_ogv,x));

eqn0 = simplify(dydx_ogv==dydx_circ);
% expr1 = simplify(subs(eqn0,{L0,Rb,Rn,rho_ogv,x0},{0.5,45/1000,12.5/1000,2.8,s1.x0(1)}));
% x_1 = vpasolve(expr1,x,[0,0.1]);

% expr2 = simplify(subs(eqn0,{L0,Rb,Rn,rho_ogv,x0},{0.5,45/1000,12.5/1000,2.8,s1.x0(2)}));



x_2 = vpasolve(expr2,x,[0,0.1]);

x0_guess = double(x_2);
rho_ogv = 2.8;
L0      = 0.5;
Rn      = 12.5/1000;
Rb      = 45/1000;
y_ogive = @(x) sqrt(rho_ogv^2 - (L0-x).^2) + Rb - rho_ogv;
% dydx_ogive = @(x) (L0-x)./sqrt(rho_ogv^2 - (L0 - x).^2);
y_circ    = @(x,x0) sqrt( Rn^2 - (x-x0).^2 );
hold on
fplot(y_ogive,[0,0.1],'r')
fplot(@(x)y_circ(x,x0_guess),[x0_guess-Rn,x0_guess+Rn],'b')
axis equal
hold off


% 
% y_circ    = @(x,x0) sqrt( Rn^2 - (x-x0).^2 );
% dydx_circ = @(x,x0) -(x-x0)./sqrt(Rn^2- (x-x0).^2);
% 
% % ogive_y = @(x,d,L) sqrt( (d*((L/d)^2+1/4))^2 -(x-L0).^2 ) - (d*((L/d)^2-1/4));
% 
% x0_guess = 0.075;
% hold on
% fplot(y_ogive,[0,0.1],'r')
% fplot(@(x)y_circ(x,x0_guess),[x0_guess-Rn,x0_guess+Rn],'b')
% axis equal