%% Ringleb symbolic stuff (v2) 09/28/2023
clc; clear; close all

syms q k a r gamma x0 y0 q0 q1 real

a_fun(q) = sqrt(1-(gamma-1)/2*q.^2);
rho_fun(q) = a_fun(q).^(2/(gamma-1));

int(1./(rho_fun(q).*q.^3),q,q0,q1)