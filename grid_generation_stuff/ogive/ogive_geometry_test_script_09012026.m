%% Example Script for ogive_geometry class (09/01/2026)
clc; clear; close all;
load("xy_ogive.mat")
x = xy_ogive(:,1);
y = xy_ogive(:,2);


%% nominal values
L0 = 500/1000;
Rb = 45.0/1000;
Rn = 12.5/1000;


%% values derived from DNS point set (uncomment to use these)

%% length from DNS points
% L0 = x(end);

%% Rb from DNS points
% Rb = y(end);

%% approximate Rn from DNS points
% p = [1, 250, 500];
% [~,~,Rn] = circ_curv3( x(p(1)), y(p(1)), ...
%                        x(p(2)), y(p(2)), ...
%                        x(p(3)), y(p(3)) );

o = ogive_geometry(L0=L0,Rn=Rn,Rb=Rb);

%% get the y coordinates from analytic relations using x-coords from DNS points
y1 = o.y_from_x(x);

% difference in the y-coordinate
e1 = abs(y1-y);

% plot
figure
hold on
plot(e1)
set(gca,'YScale','log')


%% code to generate points directly in MATLAB
% n_tip  = 298;
% n_body = 704;
% [xo,yo]  = o.get_coords(n_tip,n_body);
% hold on;
% plot(xo,yo,'r.');
% axis equal


% local function for estimating the radius of the circle from 3 points
function [kappa,c,r] = circ_curv3(x1,y1,x2,y2,x3,y3)

z1 = complex(x1,y1);
z2 = complex(x2,y2);
z3 = complex(x3,y3);

if ( ( z1 == z2 )||( z2 == z3 )||( z3 == z1 ) )
    kappa = 0;
    return
end

w = (z3 - z1)/(z2 - z1);

if ( abs(imag(w)) < eps(1) )
    kappa = 0;
    return
end

c = (z2 - z1)*(w - abs(w).^2)/(2i*imag(w)) + z1;
r = abs(z1 - c);

kappa = 1/r;

c = [real(c),imag(c)];

end