%% KT airfoil class testing (08/18/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.1;
kappa   = 0.0;
tau     = 0.0;
vinf    = 75.0;
rhoinf  = 1.0;
pinf    = 100000.0;
gamma   = 1.4;

alpha   = 10; % (degrees)
rho_ref = 1.0;
p_ref   = 100000.0;
a_ref   = sqrt(gamma*p_ref/rho_ref);

% nondimensionalize inputs
vinf   = vinf/a_ref;
rhoinf = rhoinf/rho_ref;
pinf   = pinf/(rho_ref*a_ref^2);

airfoil        = kt_airfoil(epsilon,kappa,tau);
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(alpha);

% contour coordinates
% x = [-10,0,0,-10];
% y = [-5,-5,5,5];
% src = airfoil.integrate_polygon(x,y,true);

% F = airfoil.get_flux_integrals(x(1),y(1),x(2),y(2),true) + airfoil.get_flux_integrals(x(2),y(2),x(3),y(3),true) + airfoil.get_flux_integrals(x(3),y(3),x(4),y(4),true) + airfoil.get_flux_integrals(x(4),y(4),x(1),y(1),true)
folder = fullfile(parent_dir,'MMS','kt_airfoil');


grid = read_grd_file_to_struct(fullfile(folder,'kt.grd'));
grid = grid.gblock;

folder = fullfile(parent_dir,'MMS','kt_airfoil','soln_dat_files');

DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p1_q4.dat'));
P1Q04 = DATA.ZONE.SRC_ENERGY;
DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p1_q8.dat'));
P1Q08 = DATA.ZONE.SRC_ENERGY;
DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p1_q12.dat'));
P1Q12 = DATA.ZONE.SRC_ENERGY;

DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p2_q4.dat'));
P2Q04 = DATA.ZONE.SRC_ENERGY;
DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p2_q8.dat'));
P2Q08 = DATA.ZONE.SRC_ENERGY;
DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p2_q12.dat'));
P2Q12 = DATA.ZONE.SRC_ENERGY;

DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p4_q4.dat'));
P4Q04 = DATA.ZONE.SRC_ENERGY;
DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p4_q8.dat'));
P4Q08 = DATA.ZONE.SRC_ENERGY;
DATA = tec2mat_structured_condensed(fullfile(folder,'kt-soln_p4_q12.dat'));
P4Q12 = DATA.ZONE.SRC_ENERGY;

source_line = zeros(grid.imax-1,grid.jmax-1,5);
source_curv = zeros(grid.imax-1,grid.jmax-1,5);
for j = 1:1%grid.jmax-1
    for i = 1:grid.imax-1
        x = [grid.x(i,j),grid.x(i+1,j),grid.x(i+1,j+1),grid.x(i,j+1)];
        y = [grid.y(i,j),grid.y(i+1,j),grid.y(i+1,j+1),grid.y(i,j+1)];
        source_tmp = airfoil.calculate_source(x,y,true,false);
        source_line(i,j,:) = source_tmp;
        source_curv(i,j,:) = source_tmp;
        if j == 1 % use analytic geometry
            source_tmp_curv = airfoil.calculate_source(x,y,true,true);
            source_curv(i,j,:) = source_tmp_curv;
        end
    end
end

%%
clf
hold on
plot(P1Q04(:,1)-source_curv(:,1,5),'r')
plot(P1Q08(:,1)-source_curv(:,1,5),'r--')
plot(P1Q12(:,1)-source_curv(:,1,5),'r:')

plot(P2Q04(:,1)-source_curv(:,1,5),'b')
plot(P2Q08(:,1)-source_curv(:,1,5),'b--')
plot(P2Q12(:,1)-source_curv(:,1,5),'b:')

plot(P4Q04(:,1)-source_curv(:,1,5),'g')
plot(P4Q08(:,1)-source_curv(:,1,5),'g--')
plot(P4Q12(:,1)-source_curv(:,1,5),'g:')

legend({'P1, Q4',...
        'P1, Q8',...
        'P1, Q12',...
        'P2, Q4',...
        'P2, Q8',...
        'P2, Q12',...
        'P4, Q4',...
        'P4, Q8',...
        'P4, Q12'},"NumColumns",3)