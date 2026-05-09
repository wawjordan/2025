%% Curved Wall Distance Calculation
% modifications: 04/14/2026
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
% npts = 4;
% syms xi eta
% syms xp yp zp
% syms x_(xi,eta)
% syms y_(xi,eta)
% syms z_(xi,eta)
% syms L_(xi,eta)
% syms num_(xi,eta)
% 
% 
% tmp = diff(num_./L_,xi);
% % L = ( (x_-xp).^2 + (y_-yp).^2 + (z_-zp).^2 ).^(1/2);
% % 
% % f1 = diff(L,xi);
% % f2 = diff(L,eta);
% 
% % a = diff(f1,xi);
% % b = diff(f1,eta);
% % c = diff(f2,xi);
% % d = diff(f2,eta);
% % det = a.*d - b.*c;

syms A A1 A2
syms B B1 B2
syms L

detJ = A1*B2*L^4 - A1*L^2*B^2 - B2*L^2*A^2 + A^2*B^2 -( A2*B1*L^4 - A2*L^2*A*B - B1*L^2*A*B + A^2*B^2);
detJ = simplify(subs(detJ,B1,A2));

Jinv11 = B2*L^2 - B^2;
Jinv12 = A*B - A2*L^2;
Jinv22 = A1*L^2 - A^2;

up1 = (Jinv11*(A/L) + Jinv12*(B/L) )/detJ;
up2 = (Jinv12*(A/L) + Jinv22*(B/L) )/detJ;

pretty(up1)
pretty(simplify(up1))