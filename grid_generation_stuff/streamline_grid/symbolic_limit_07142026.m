%% Joukowski Airfoil Limits
% 07/14/2026
clc; clear; close all;
syms l mu real
syms epsilon kappa alpha
syms zeta complex
% z = joukowski_map(zeta,l);
% w = airfoil_velocity(x,y,l,epsilon,kappa,tau,alpha,vinf);

function w = w_from_zeta(zeta,l,mu,c1,c2,c3)
% a    = l*sqrt( (1+epsilon)^2 + kappa^2 );
% beta = asin(l*kappa/a);
% mu   = l*(-epsilon + kappa*1i);
% c1 = exp(-1i*alpha)
% c2 = 2*1i*a*sin(alpha+beta)
% c3 = - a^2*exp(1i*alpha)
dz = dmap(zeta,l);
% wtilde = exp(-1i*alpha) + 2*1i*a*sin(alpha+beta)./(zeta-mu) - a^2*exp(1i*alpha)./(zeta-mu).^2;
wtilde = c1 + c2./(zeta-mu) +c3./(zeta-mu).^2;
w = simplify( wtilde./dz );
end

function z = map(zeta,l)
z = zeta + l^2./zeta;
end

function z = dmap(zeta,l)
z = 1 - l^2./zeta.^2;
end