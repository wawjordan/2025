%% Ringleb Flow Grid Plotting
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\ringleb_flow\ring\';
% subfolder = '_ring0257x1025\';
subfolder = '_ring0033x0129\';
file      = 'ring.grd';


GRID = read_grd_file( fullfile(folder,subfolder,file) );
GRID = GRID.gblock;
GRID.x = GRID.x(:,:,1);
GRID.y = GRID.y(:,:,1);
hold on;
plot(GRID.x,GRID.y,'k')
plot(GRID.x.',GRID.y.','k')
axis equal

% X1 = GRID.x;
% Y1 = GRID.y;
% max_iter = 1;
% tol = 1e-1;
% offset = floor(GRID.jmax/2)-1;

% [X2,Y2,R] = winslow_eqns_fd(X1,Y1,max_iter,tol);
% [X2,Y2,R] = winslow_eqns_fd(X1(:,1+offset:GRID.jmax-offset),Y1(:,1+offset:GRID.jmax-offset),max_iter,tol);
% X1(:,1+offset:GRID.jmax-offset) = X2;
% Y1(:,1+offset:GRID.jmax-offset) = Y2;
% plot(X1,Y1,'r')
% plot(X1.',Y1.','r')


% clf;
% normal = [0.2850695464146078,0.95850683550351889];
% s = [normal(2),-normal(1)];
% v1 = [0.38656862376259332,-0.12855270379788061];
% v2 = [0.39399166546473552,-0.10359375500374032];
% hold on;
% quiver(0,0,s(1),s(2),'k')
% quiver(0,0,normal(1),normal(2),'g')
% quiver(0,0,v1(1),v1(2),'r')
% quiver(0,0,v2(1),v2(2),'b')
% axis equal