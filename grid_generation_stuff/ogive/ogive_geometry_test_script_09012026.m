%% Example Script for ogive_geometry class (09/01/2026)
clc; clear; close all;

o = ogive_geometry(L0=0.5,Rn=12.5/1000,Rb=45.0/1000);
n_tip  = 298;
n_body = 704;
[x,y]  = o.get_coords(n_tip,n_body);

hold on;
plot(x,y,'r.');
axis equal

%% uncomment to plot the DNS points as well
% load("xy_ogive.mat")
% plot(xy_ogive(:,1),xy_ogive(:,2),'k.');
% plot(x,y,'ro');