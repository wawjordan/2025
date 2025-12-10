%% KT airfoil without complex arithmetic
clc; clear; close all;

l          = 1;
epsilon    = 0.1;
kappa      = -0.0;
% tau        = deg2rad(0.174532925199428);
tau        = 0.1;

% theta = linspace(0,3*pi/2,101);
theta = linspace(-1,0.5,101);

[x,y] = kt_from_theta(l,epsilon,kappa,tau,theta);

load("kt_x.mat")
load("kt_y.mat")
plot(xtmp,ytmp,'k--')
hold on;
plot(x,y,'r')
axis equal


function [x,y] = kt_from_theta(l,epsilon,kappa,tau,t)
% calculate the non-dimensional radial coordinate in the circle plane
% (R = r/b), from eq. 5 (sign error?)
b = l;
m = 2 - tau/pi;

% theta = 2*pi - theta;
theta = pi*(-t);

F = epsilon/b;
G = kappa/b;
R = F*cos(theta) + G*sin(theta) + sqrt( 1 + 2*F + ( F*cos(theta) + G*sin(theta) ).^2 );

% calculate the corresponding geometric segments R1 and R2 and their ratio,
% t, from eq. 4
R1 = max(sqrt( 1 + R.^2 + 2*R.*cos(theta) ),1e-14);
R2 = max(sqrt( 1 + R.^2 - 2*R.*cos(theta) ),1e-14);
t = R2./R1;

% now calculate (phi_2 - phi_1) (also eq. 4)
% d_phi = atan( 2*sin(theta) / ( R - (1/R) ) );
% if (R<1)
%     d_phi = d_phi - pi;
% end
d_phi = atan2( 2*sin(theta), R - (1./R) );

% calculate the divisor for defining xi and eta (from eq. 3)
% note that the zeta (xi,i*eta) <=> z (x,i*y) notation is reversed here
% (x,y) are (xi,eta) which are the physical airfoil coordinates
den = 1 + t.^(2*m) - 2*(t.^m).*cos(m*d_phi);

x = m*b*( ( 1 - t.^(2*m) )./den );
x = -x;
y = m*b*( ( 2*(t.^m).*sin(m*d_phi) )./den );

% xLE = m*b*( (1+F).^(2*m) - F.^(2*m) )./( ( (1+F).^m - F.^m).^2 );
% x = xLE - x;
end