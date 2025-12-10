%% Ringleb symbolic stuff
clc; clear; close all

syms q k a r gamma x0 y0

syms c(q) rho(q) J(q)

a_fun(q) = sqrt(1-(gamma-1)/2*q.^2);
rho_fun(a) = a.^(2/(gamma-1));
% p_fun(a) = 1/gamma * a.^(2*gamma/(gamma-1));
J_fun(a) = 1./a + 1./(3*a.^3) + 1./(5*a.^5) - (1/2)*log( (1+a)./(1-a) );

dadq(q) = simplify(diff(a_fun(q),q));
dadq(q) = subs(dadq(q),a_fun(q),a);

drhoda(a) = simplify(diff(rho_fun(a),a));

dJda(a) = simplify(diff(J_fun(a),a));

x_fun(q,k) = 1./(2*rho).*( 2./k.^2 - 1./q.^2 ) - J/2;
y_fun(q,k) = sqrt( 1 - (q./k).^2 )./(k.*rho.*q);

dxdq = diff(x_fun(q,k),q);
dxdq = subs(dxdq, diff(J(q), q) , dJda(a)*dadq(q) );
dxdq = subs(dxdq,diff(rho(q), q), drhoda(a)*dadq(q));
dxdq = subs(dxdq,rho(q), rho_fun(a));
dxdq = subs(dxdq,a,a_fun(q));
dxdq = simplify(dxdq);

dxdk = diff(x_fun(q,k),k);
% dxdk = subs(dxdk,rho(q), rho_fun(a));
% dxdk = subs(dxdk,a,a_fun(q));
% dxdk = simplify(dxdk);


% dxdq = matlabFunction(diff(x,q));
% dxdk = matlabFunction(diff(x,k));
% dydq = matlabFunction(diff(y,q));
% dydk = matlabFunction(diff(x,k));
