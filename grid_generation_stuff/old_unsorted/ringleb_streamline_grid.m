function [x,y,q,k] = ringleb_streamline_grid(kmin,kmax,q0,Nq,Nk,spacing)
k = spacing{1}(kmin,kmax,Nk);
q = zeros(Nq,Nk);
for j = 1:Nk
    % decide spacing of potential lines
    theta_min = asin(q0/k(j));
    theta_max = pi/2;
    phi_start = k(j)/q0  *cos(theta_min);
    phi_end   = k(j)/k(j)*cos(theta_max);
    phi = spacing{2}(phi_start,phi_end,Nq);
    for i = 1:Nq
        q(i,j) = sqrt(k(j)^2/(1+phi(i)^2));
    end
end
k = repmat(k(:).',[Nq,1]);
x = x_fun(q,k);
y = y_fun(q,k);
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

% function val = kmin;    val = 1.0;                         end
% function val = kmax;    val = 1.5;                         end
% function val =   q0;    val = 0.8;                         end

% function val = kmin;    val = 1.1;                         end
% function val = kmax;    val = 1.6;                         end
% function val =   q0;    val = 1.05;                         end

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