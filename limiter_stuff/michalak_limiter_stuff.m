clear; close all; clc;

syms r

% min(1,r) that we want to approximate
phi0 = piecewise( r < 0, 0, (0 <= r) & (r <= 1.0),r,r > 1.0, 1);
dphi0_dr = diff(phi0,r);

phi02 = piecewise( 2*r < 0, 0, (0 <= 2*r) & (2*r <= 1.0),2*r,2*r > 1.0, 1);
dphi02_dr = diff(phi02,r);

% Michalak's polynomial function for yt=1.5
% phi1 = piecewise( r < 0, 0, (0 <= r) & (r <= 1.5),r-(4/27)*r.^3,r > 1.5, 1);
phi1 = piecewise( 2*r < 0, 0, (0 <= 2*r) & (2*r <= 1.5),(2*r)-(4/27)*(2*r).^3,2*r > 1.5, 1);
dphi1_dr = diff(phi1,r);

% proposed modification
phi2 = piecewise( r < 0, 0, (0 <= r) & (r <= 1.0),(3/2)*r-(1/2)*r.^3,r > 1.0, 1);
dphi2_dr = diff(phi2,r);

figure(1);
hold on;
fplot(phi0,[-1,2],'k')
fplot(phi02,[-1,2],'k--')
fplot(phi1,[-1,2],'r')
fplot(phi2,[-1,2],'b')
legend('min(1,r)','min(1,2*r)','Michalak','Proposed')

% figure(2);
% hold on;
% fplot(dphi02_dr,[-1,2],'k--')
% fplot(dphi1_dr,[-1,2],'r')
% fplot(dphi2_dr,[-1,2],'b')