%% Example Script for ogive_geometry class (09/01/2026)
clc; clear; close all;

o = ogive_geometry(L0=0.5,Rn=12.5/1000,Rb=45.0/1000);
% n_tip  = 298;
% n_body = 704;
% [x,y]  = o.get_coords(n_tip,n_body);
% 
% hold on;
% plot(x,y,'r.');
% axis equal

%% uncomment to plot the DNS points as well
load("xy_ogive.mat")
x = xy_ogive(:,1);
y = xy_ogive(:,2);
mask = x<o.xt;

x = x(mask);
y = y(mask);


N = numel(x);

p = [1, 250, 500];

[~,c,r] = circ_curv3( x(p(1)), y(p(1)), ...
                      x(p(2)), y(p(2)), ...
                      x(p(3)), y(p(3)) );


% plot(xy_ogive(:,1),xy_ogive(:,2),'k.');
% plot(x,y,'ro');

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