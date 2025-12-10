function [x,y,q,k] = ringleb_transfinite_spline_grid(kmin,kmax,q0,Nq,Nk,spacing)
F = construct_transfinite_spline_map(kmin,kmax,q0,Nq,Nk,spacing);
[qt,kt] = ndgrid(linspace(-1,1,Nq),linspace(-1,1,Nk));
q = qt;
k = kt;
x = F.x(qt,kt);
y = F.y(qt,kt);
end

function a = sound_speed(q)
a = sqrt(one - half*gm1*q.^2);
end

function d = dadq_fun(q)
d = -half*gm1*q./sound_speed(q);
end

function rho = density(q)
a = sound_speed(q);
rho = a.^(two*xgm1);
end

function d = drhodq_fun(q)
% d = -q.*(one - half*gm1*q.^2).^(xgm1-one);
a = sound_speed(q);
drhoda = (two*a.^(two*xgm1 - one))*xgm1;
d = drhoda.*dadq_fun(q);
end

function J = Jquantity(q)
a = sound_speed(q);
J = one./a + one./(three*a.^3) + one./(five*a.^5) - half*log((one+a)./(one-a));
end

function d = dJda_fun(a)
d = one./(a.^6.*(a.^2-one));
end

function x = x_fun(q,k)
rho = density(q);
J   = Jquantity(q);
x = half./rho .* (two./k.^2 - one./q.^2) - half*J;
end

function d = dxdq_fun(q,k)
rho = density(q);
drhodq = drhodq_fun(q);
a = sound_speed(q);
dadq = dadq_fun(q);
dJda = dJda_fun(a);

% d1 = -one./(k.^2 .* rho.^2) .* drhodq;
% d2 = half*( one./(q.^2 .* rho.^2) .* drhodq + two./(rho.*q.^3));
% d3 = -half.*dJda.*dadq;
% d = d1 + d2 + d3;
d1 = one./(rho.*q.^3);
d2 = -half*drhodq.*(two./k.^2 - one./q.^2);
d3 = -half.*dJda.*dadq;
d = d1 + d2 + d3;
end

function d = dxdk_fun(q,k)
rho = density(q);
d = -two./ (rho .* k.^3);
end

function y = y_fun(q,k)
rho = density(q);
y = sqrt(one - (q./k).^2) ./ (k.*rho.*q);
end

function d = dydq_fun(q,k)
rho = density(q);
drhodq = drhodq_fun(q);
d1 = sqrt(one - (q./k).^2);
d2 = -one./( (k.^3 .* rho) .* d1 );
d3 = -( d1./(k .* q    .* rho.^2) ).*drhodq;
d4 = -( d1./(k .* q.^2 .* rho   ) );
d = d2 + d3 + d4;
end

function d = dydk_fun(q,k)
rho = density(q);
d1 = sqrt(one - (q./k).^2);
d2 = q ./ ( k.^4 .* rho .* d1 );
d3 = - d1./(k.^2 .* q .* rho );
d = d2 + d3;
end




function [x,y] = inflow(q0,k)
q = q0;
k2 = reparameterize_potential_line(q,k);
x = x_fun(q,k2);
y = y_fun(q,k2);
end

function [x,y] = inner_wall(q,kmax)
k = kmax;
q2 = reparameterize_streamline(q,k);
x = x_fun(q2,k);
y = y_fun(q2,k);
end

function [x,y] = outer_wall(q,kmin)
k = kmin;
q2 = reparameterize_streamline(q,k);
x = x_fun(q2,k);
y = y_fun(q2,k);
end

function F = construct_transfinite_spline_map(kmin,kmax,q0,Nq,Nk,spacing)




% inner wall - q = [q0,kmax], k = [kmax,kmax]
XY1 = zeros(Nq,2);
q = spacing{1}(q0,kmax,Nq);
[XY1(:,1),XY1(:,2)] = inner_wall(q,kmax);

% symmetry - q = [kmin,kmax], k = [kmin,kmax]
XY2 = zeros(Nk,2);
% XY2(:,1) = spacing{2}(x_fun(kmin,kmin),x_fun(kmax,kmax),Nk);
% XY2(:,2) = spacing{2}(y_fun(kmin,kmin),y_fun(kmax,kmax),Nk);
% Note this:
XY2(:,1) = spacing{2}(x_fun(kmax,kmax),x_fun(kmin,kmin),Nk);
% XY2(:,2) = spacing{2}(y_fun(kmax,kmax),y_fun(kmin,kmin),Nk);
XY2(:,2) = zeros(Nk,1);

% outer wall - q = [q0,kmin], k = [kmin,kmin]
XY3 = zeros(Nq,2);
q = spacing{1}(q0,kmin,Nq);
[XY3(:,1),XY3(:,2)] = outer_wall(q,kmin);

% inflow     - q = [q0,q0], k = [kmin,kmax]
XY4 = zeros(Nk,2);
k = spacing{2}(kmax,kmin,Nk);
[XY4(:,1),XY4(:,2)] = inflow(q0,k);

F = transfinite_spline_interpolant_2D(XY1,XY2,XY3,XY4);

end

function q2 = reparameterize_streamline(q1,k)
N = length(q1);

darcL = @(q) sqrt(dxdq_fun(q,k).^2 + dydq_fun(q,k).^2);
% arc length approximation
L = integral( @(q) darcL(q), q1(1), q1(N),'ArrayValued',true );

% arc length parameterization
L1 = q1(N) - q1(1);
q2 = q1;
for j = 2:N-1
    f2 = @(t) q1(j) - q1(j-1) - (L1/L)*integral( @(z) darcL(z), q2(j-1), t,'ArrayValued',true );
    a1 = q1(1); % needs to be a larger value than the previous one
    b1 = q1(N);  %
    % x1 = linspace( a1, b1, 10 );
    % hold on
    % plot(x1,arrayfun(f2,x1))
    q2(j) = fzero( f2, [a1,b1]);
end

end

function k2 = reparameterize_potential_line(q,k1)
N = length(k1);

darcL = @(k) sqrt(dxdk_fun(q,k).^2 + dydk_fun(q,k).^2);
% arc length approximation
L = integral( @(k) darcL(k), k1(1), k1(N),'ArrayValued',true );

% arc length parameterization
L1 = k1(N) - k1(1);
k2 = k1;
for j = 2:N-1
    f2 = @(t) k1(j) - k1(j-1) - (L1/L)*integral( @(z) darcL(z), k2(j-1), t,'ArrayValued',true );
    a1 = k1(1);
    b1 = k1(N);  %
    % x1 = linspace( a1, b1, 10 );
    % hold on
    % plot(x1,arrayfun(f2,x1))
    k2(j) = fzero( f2, [a1,b1]);
end

end

% function val = kmin;    val = 0.7;                         end
% function val = kmax;    val = 1.5;                         end
% function val =   q0;    val = 0.4;                         end
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