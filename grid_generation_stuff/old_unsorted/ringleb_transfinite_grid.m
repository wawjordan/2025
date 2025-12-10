function [x,y,q,k] = ringleb_transfinite_grid(kmin,kmax,q0,Nq,Nk,spacing)
F = construct_transfinite_map(kmin,kmax,q0);
[qt,kt] = ndgrid(spacing{1}(-1,1,Nq),spacing{2}(-1,1,Nk));
q = qt;
k = kt;
x = F.x(qt,kt);
y = F.y(qt,kt);
end

function q = qmap1(x,~,kmax,q0)
q = 0.5*( (1-x)*q0  + (1+x)*kmax );
end

function q = qmap2(x,kmin,~,q0)
q = 0.5*( (1-x)*q0 + (1+x)*kmin );
end

function k = kmap(x,kmin,kmax,~)
k = 0.5*( (1-x)*kmin + (1+x)*kmax );
end

function a = sound_speed(q)
a = sqrt(one - half*gm1*q.^2);
end

function rho = density(q)
a = sound_speed(q);
rho = a.^(two*xgm1);
end

function J = Jquantity(q)
a = sound_speed(q);
J = one./a + one./(three*a.^3) + one./(five*a.^5) - half*log((one+a)./(one-a));
end

function x = x_fun(q,k)
rho = density(q);
J   = Jquantity(q);
x = half./rho .* (two./k.^2 - one./q.^2) - half*J;
end

function y = y_fun(q,k)
rho = density(q);
y = sqrt(one - (q./k).^2) ./ (k.*rho.*q);
end

function F = construct_transfinite_map(kmin,kmax,q0)

% inner wall - q = [q0,kmax], k = [kmax,kmax]
G01.x = @(qt) qmap1(qt,kmin,kmax,q0);
G01.y = @(qt) 0*qt + kmap(1,kmin,kmax,q0);

% symmetry - q = [kmin,kmax], k = [kmin,kmax]
G02.x = @(kt) kmap(kt,kmin,kmax,q0);
G02.y = @(kt) kmap(kt,kmin,kmax,q0);

% outer wall - q = [q0,kmin], k = [kmin,kmin]
G03.x = @(qt) qmap2(qt,kmin,kmax,q0);
G03.y = @(qt) 0*qt + kmap(-1,kmin,kmax,q0);

% inflow     - q = [q0,q0], k = [kmin,kmax]
G04.x = @(kt) 0*kt + qmap2(-1,kmin,kmax,q0);
G04.y = @(kt) kmap(kt,kmin,kmax,q0);

F1 = transfinite_interpolant(G01,G02,G03,G04);

% outer wall - q = [q0,kmin], k = [kmin,kmin]
G1.x = @(qt) x_fun(F1.x(qt,-1),F1.y(qt,-1));
G1.y = @(qt) y_fun(F1.x(qt,-1),F1.y(qt,-1));

% symmetry - q = [kmin,kmax], k = [kmin,kmax]
G2.x = @(kt) x_fun(F1.x(1,kt),F1.y(1,kt));
G2.y = @(kt) y_fun(F1.x(1,kt),F1.y(1,kt));

% inner wall - q = [q0,kmax], k = [kmax,kmax]
G3.x = @(qt) x_fun(F1.x(qt,1),F1.y(qt,1));
G3.y = @(qt) y_fun(F1.x(qt,1),F1.y(qt,1));

% inflow     - q = [q0,q0], k = [kmin,kmax]
G4.x = @(kt) x_fun(F1.x(-1,kt),F1.y(-1,kt));
G4.y = @(kt) y_fun(F1.x(-1,kt),F1.y(-1,kt));

F = transfinite_interpolant(G1,G2,G3,G4);
end

% function val = zero;    val = 0.0;                         end
function val = one;     val = 1.0;                         end
function val = two;     val = 2.0;                         end
function val = three;   val = 3.0;                         end
% function val = four;    val = 4.0;                         end
function val = five;    val = 5.0;                         end
function val = half;    val = one/two;                     end

function val = gamma;   val = 1.4;                         end
function val = gm1;     val = gamma - one;                 end
% function val = gp1;     val = gamma + one;                 end
% function val = xg;      val = one / gamma;                 end
function val = xgm1;    val = one / (gamma - one);         end
% function val = xgp1;    val = one / (gamma + one);         end
% function val = xg2m1;   val = one / (gamma*gamma - one);   end
% function val = gxgm1;   val = gamma / (gamma - one);       end
% function val = gxg2m1;  val = gamma / (gamma*gamma - one); end
% function val = gm1xgp1; val = gm1/gp1;                     end
% function val = gp1xgm1; val = gp1/gm1;                     end