%% Ringleb Flow prototyping
clc; clear; close all;

zero; one; two; three; five; half;
gamma; gm1; gp1; xg; xgm1; xgp1; xg2m1; gxgm1; gxg2m1; gm1xgp1; gp1xgm1;
kmin; kmax; q0;


spacing1 = {@linspace,@linspace};
% spacing1 = {@(xmin,xmax,N)mygeomspace( xmin, xmax, 1.00, N ),@(xmin,xmax,N)mygeomspace( xmin, xmax, 0.99, N )};
% Nq = 129;
% Nk = 65;
Nq = 17;
Nk = 17;

% [xt,yt,qt,kt] = transfinite_grid(kmin,kmax,q0,Nq,Nk,spacing1);
[xt,yt,qt,kt] = transfinite_spline_grid(kmin,kmax,q0,Nq,Nk,spacing1);
[xs,ys,qs,ks] = streamline_grid(   Nq,Nk,spacing1);
%%
hold on; axis equal
plot(xs,ys,'k'); plot(xs.',ys.','k');
plot(xt,yt,'r'); plot(xt.',yt.','r');
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','FunctionTolerance',1e-10,'MaxFunctionEvaluations',10000,'StepTolerance',1e-17,'Display','none');
% [q,k] = qk_map(qguess,kguess,kmin,kmax,q0,x,y,tol,maxiter)
tol = 1e-10;
maxiter = 1000;
res = 0*qs;
niters = 0*qs;
xst = 0*xs;
yst = 0*ys;

for j = 1:Nk
    for i = 1:Nq
        % qguess = (kmin + q0)/2;
        % kguess = (kmax + kmin)/2;
        qguess = qs(i,j);
        kguess = ks(i,j);
        x_set = xt(i,j);
        y_set = yt(i,j);
        [qtmp,ktmp,restmp,niter] = qk_map(qguess,kguess,kmin,kmax,q0,x_set,y_set,tol,maxiter);
        qt(i,j) = qtmp;
        kt(i,j) = ktmp;
        xst(i,j) = x_fun(qtmp,ktmp);
        yst(i,j) = y_fun(qtmp,ktmp);
        res(i,j) = restmp;
        niters(i,j) = niter;
    end
end
plot(xst,yst,'b'); plot(xst.',yst.','b');


% qt2 = 0*qt;
% kt2 = 0*kt;
% xst2 = xst*0;
% yst2 = yst*0;
% dxy = @(x,xq,yq) (x_fun(x(1),x(2)) - xq).^2 + (y_fun(x(1),x(2))-yq).^2;
% 
% for j = 1:Nk
%     for i = 1:Nq
%         qguess = (kmin + q0)/2;
%         kguess = (kmax + kmin)/2;
%         x_set = xt(i,j);
%         y_set = yt(i,j);
%         tmp = fsolve(@(x)dxy(x,x_set,y_set),[qguess,kguess],options);
%         qtmp = tmp(1);
%         ktmp = tmp(2);
%         qt2(i,j) = qtmp;
%         kt2(i,j) = ktmp;
%         xst2(i,j) = x_fun(qtmp,ktmp);
%         yst2(i,j) = y_fun(qtmp,ktmp);
%     end
% end
% 
% plot(xst2,yst2,'g'); plot(xst2.',yst2.','g');
% xq = [1,1];
% dxy = @(x) (x_fun(x(1),x(2)) - xq(1)).^2 + (y_fun(x(1),x(2))-xq(2)).^2;
% 
% qk = fsolve(dxy,[q0,kmax]);
% dxy(qk)
% 
% plot(x_fun(qk(1),qk(2)),y_fun(qk(1),qk(2)),'k.')
% labels ={length(x),1};
% for i = 1:length(x)
%     labels{i} = sprintf('%d',i);
% end
% text(x,y,labels)

function q = qmap1(x)
q = 0.5*( (1-x)*q0 + (1+x)*kmax );
end

function q = qmap2(x)
q = 0.5*( (1-x)*q0 + (1+x)*kmin );
end

function k = kmap(x)
k = 0.5*( (1-x)*kmin + (1+x)*kmax );
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
den = a.^6.*(a.^2-one);
den = den*(den~=0) + eps(one)*(den==0);
d = one./(den);
end

function theta = theta_fun(q,k)
theta = asin(q./k);
end

function psi = stream_fun(q,k)
theta = theta_fun(q,k);
psi = (one./q).*sin(theta);
end

function u = uvel(q,k)
theta = theta_fun(q,k);
u = q.*cos(theta);
end

function v = vvel(q,k)
theta = theta_fun(q,k);
v = q.*sin(theta);
end

function p = pressure(q)
a = sound_speed(q);
p = xg*a.^(two*gxgm1);
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
arg = one - (q./k).^2;
arg = arg*(arg>0) + eps(one)*(arg<=0);
d1 = sqrt(arg);
d2 = -one./( (k.^3 .* rho) .* d1 );
d3 = -( d1./(k .* q    .* rho.^2) ).*drhodq;
d4 = -( d1./(k .* q.^2 .* rho   ) );
d = d2 + d3 + d4;
end

function d = dydk_fun(q,k)
rho = density(q);
arg = one - (q./k).^2;
arg = arg*(arg>0) + eps(one)*(arg<=0);
d1 = sqrt(arg);
d2 = q ./ ( k.^4 .* rho .* d1 );
d3 = - d1./(k.^2 .* q .* rho );
d = d2 + d3;
end




function [x,y] = inflow(k)
q = q0;
k2 = reparameterize_potential_line(q,k);
x = x_fun(q,k2);
y = y_fun(q,k2);
end

function [x,y] = inner_wall(q)
k = kmax;
q2 = reparameterize_streamline(q,k);
x = x_fun(q2,k);
y = y_fun(q2,k);
end

function [x,y] = outer_wall(q)
k = kmin;
q2 = reparameterize_streamline(q,k);
x = x_fun(q2,k);
y = y_fun(q2,k);
end

function F = construct_transfinite_spline_map(kmin,kmax,q0,Nq,Nk,spacing)




% inner wall - q = [q0,kmax], k = [kmax,kmax]
XY1 = zeros(Nq,2);
q = spacing{1}(q0,kmax,Nq);
[XY1(:,1),XY1(:,2)] = inner_wall(q);

% symmetry - q = [kmin,kmax], k = [kmin,kmax]
XY2 = zeros(Nk,2);
% XY2(:,1) = spacing{2}(x_fun(kmin,kmin),x_fun(kmax,kmax),Nk);
% XY2(:,2) = spacing{2}(y_fun(kmin,kmin),y_fun(kmax,kmax),Nk);
% Note this:
XY2(:,1) = spacing{2}(x_fun(kmax,kmax),x_fun(kmin,kmin),Nk);
XY2(:,2) = spacing{2}(y_fun(kmax,kmax),y_fun(kmin,kmin),Nk);

% outer wall - q = [q0,kmin], k = [kmin,kmin]
XY3 = zeros(Nq,2);
q = spacing{1}(q0,kmin,Nq);
[XY3(:,1),XY3(:,2)] = outer_wall(q);

% inflow     - q = [q0,q0], k = [kmin,kmax]
XY4 = zeros(Nk,2);
k = spacing{2}(kmax,kmin,Nk);
[XY4(:,1),XY4(:,2)] = inflow(k);

F = transfinite_spline_interpolant_2D(XY1,XY2,XY3,XY4);

end

function F = construct_transfinite_map()

% inner wall - q = [q0,kmax], k = [kmax,kmax]
G01.x = @(qt) qmap1(qt);
G01.y = @(qt) 0*qt + kmap(1);

% symmetry - q = [kmin,kmax], k = [kmin,kmax]
G02.x = @(kt) kmap(kt);
G02.y = @(kt) kmap(kt);

% outer wall - q = [q0,kmin], k = [kmin,kmin]
G03.x = @(qt) qmap2(qt);
G03.y = @(qt) 0*qt + kmap(-1);

% inflow     - q = [q0,q0], k = [kmin,kmax]
G04.x = @(kt) 0*kt + qmap2(-1);
G04.y = @(kt) kmap(kt);

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

function [x,y,q,k] = transfinite_grid(~,~,~,Nq,Nk,spacing)
F = construct_transfinite_map();
[qt,kt] = ndgrid(spacing{1}(-1,1,Nq),spacing{2}(-1,1,Nk));
q = qt;
k = kt;
x = F.x(qt,kt);
y = F.y(qt,kt);
end

function [x,y,q,k] = transfinite_spline_grid(kmin,kmax,q0,Nq,Nk,spacing)
F = construct_transfinite_spline_map(kmin,kmax,q0,Nq,Nk,spacing);
[qt,kt] = ndgrid(linspace(-1,1,Nq),linspace(-1,1,Nk));
q = qt;
k = kt;
x = F.x(qt,kt);
y = F.y(qt,kt);
end

function [x,y] = stream_line(q)
k = kmax;
x = x_fun(q,k);
y = y_fun(q,k);
end


function [x,y,q,k] = streamline_grid(Nq,Nk,spacing)
k = spacing{1}(kmin,kmax,Nk);
q = zeros(Nq,Nk);
for j = 1:Nk
    % q1 = spacing{2}(q0,k(j),Nq);
    % q2 = reparameterize_streamline(q1,k(j));
    % q2 = reparameterize_streamline_2(q1,k(j));
    % q(:,j) = q2;

    % decide spacing of potential lines
    theta_min = asin(q0/k(j));
    theta_max = pi/2;
    phi_start = k(j)/q0  *cos(theta_min);
    phi_end   = k(j)/k(j)*cos(theta_max);
    phi = spacing{2}(phi_start,phi_end,Nq);
    % fun = @(q,p) p - (k(j)./q).*cos(asin(q./k(j)));
    for i = 1:Nq
        % q(i,j) = fzero(@(q)fun(q,phi(i)),[q0,k(j)]);
        % q1 = fzero(@(q)fun(q,phi(i)),[q0,k(j)]);
        % q2 = sqrt(k(j)^2/(1+phi(i)^2));
        q(i,j) = sqrt(k(j)^2/(1+phi(i)^2));
    end
end
k = repmat(k(:).',[Nq,1]);
x = x_fun(q,k);
y = y_fun(q,k);
end

% function [x,y] = boundaries(Nq,Nk,spacing)
% 
% k = spacing{1}(kmin,kmax,Nk).';
% q = spacing{2}(q0,kmin,Nq).';
% q4 = spacing{2}(q0,kmax,Nq);
% 
% [x2,y2] = outer_wall(q2);
% [x3,y3] = inflow(k3);
% [x4,y4] = inner_wall(q4);
% 
% x = [flip(x2(:));x3(:);x4(:)];
% y = [flip(y2(:));y3(:);y4(:)];
% end

% function x = cgb_space(x_start,x_end,N)
% d = abs(x_end-x_start);
% [x0,~] = ChebGaussLob(N);
% x0 = (x0+1)/2;
% if x_start < x_end
%     x = x_start + x0*d;
% else
%     x = x_start - x0*d;
% end
% end

function x = cos_space(x_start,x_end,N)
d = x_end-x_start;
x0 = cos(linspace(-pi,0,N));
x0 = (x0+1)/2;
x = x_start + x0*d;
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

% function q2 = reparameterize_streamline_2(q1,k)
% N = length(q1);
% 
% darcL = @(q) 1+ tanh(q+5);
% q2 = weight_grid(q1(1),q1(N),N,darcL);
% 
% end

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


% function [qnp1,knp1] = qk_update(qn,kn,x,y)
% 
% dxn = x_fun(qn,kn) - x;
% dyn = y_fun(qn,kn) - y;
% 
% dxdq = dxdq_fun(q,k);
% dxdk = dxdk_fun(q,k);
% dydq = dydq_fun(q,k);
% dydk = dydk_fun(q,k);
% 
% Jinv = dxdq.*dydk - dxdk.*dydq;
% 
% qnp1 = qn - Jinv.*( dydk.*dxn - dxdk.*dyn );
% knp1 = kn - Jinv.*( dxdq.*dyn - dydq.*dxn );
% 
% end

function [q,k,R,N] = qk_map(qguess,kguess,kmin,kmax,q0,x,y,tol,maxiter)
q = qguess;
k = kguess;
% get initial deviation
dxn = x_fun(q,k) - x;
dyn = y_fun(q,k) - y;

% dt  = 10;
% dtmin = 1.0e-5;
% dtmax = 1.0e+5;
small = 1e-6;

err0 = sqrt(dxn.^2 + dyn.^2) + small;

R = one;
for i = 1:maxiter
    N = i;
    if (R<tol)
        break
    end
    dxdq = dxdq_fun(q,k);
    dxdk = dxdk_fun(q,k);
    dydq = dydq_fun(q,k);
    dydk = dydk_fun(q,k);
    Jinv = 1./(dxdq.*dydk - dxdk.*dydq);
    % Jinv = 1./((dxdq+1/dt).*(dydk+1/dt) - dxdk.*dydq);

    % tmp = [q;k] - Jinv*[1/dt + dydk, - dxdk; -dydq, 1/dt + dxdq]*[dxn;dyn];
    % tmp = [q;k] - Jinv*[dydk, - dxdk; -dydq, dxdq]*[dxn;dyn];
    % qnew = tmp(1);
    % knew = tmp(2);

    dq = Jinv.*( dydk.*dxn - dxdk.*dyn );
    dk = Jinv.*( dxdq.*dyn - dydq.*dxn );
    % dq = Jinv.*( (dydk + 1/dt).*dxn - dxdk.*dyn );
    % dk = Jinv.*( (dxdq + 1/dt).*dyn - dydq.*dxn );
    

    in_bounds = false;
    UR_factor = 1;
    while ~in_bounds
        qnew = q - UR_factor*dq;
        knew = k - UR_factor*dk;

        % tmp = [q;k] - UR_factor*[dxdq, dxdk; dydq, dydk]\[q;k];
        % qnew = tmp(1);
        % knew = tmp(2);
        in_bounds = bounds_check(qnew,knew,q0,kmin,kmax);
        UR_factor = UR_factor/2;
        if (isnan(qnew)||isnan(knew))
            error('nan in qk_map');
        end
    end
    q = qnew;
    k = knew;

    
    dxn = x_fun(q,k) - x;
    dyn = y_fun(q,k) - y;
    % plot(x+dxn,y+dyn,'o')

    R = sqrt(dxn.^2 + dyn.^2)./err0;
    % R = cond(Jinv*[dydk, - dxdk; -dydq, dxdq]);

    % dt = dt/R;
    % dt = max(dt,dtmin);
    % dt = min(dt,dtmax);
end
if (R>tol)
    warning('tolerance not met in qk_map')
end
% plot(x+dxn,y+dyn,'o')
end

function in_bounds = bounds_check(q,k,q0,kmin,kmax)
left_check   = ( q >= q0 );
right_check  = ( q <= k );
bottom_check = ( k >= kmin ) & (k >= q);
top_check    = ( k <= kmax );

in_bounds = left_check & right_check & bottom_check & top_check;
end






function val = kmin;    val = 0.7;                         end
function val = kmax;    val = 1.5;                         end
function val =   q0;    val = 0.4;                         end
% function val = kmin;    val = 1.0;                         end
% function val = kmax;    val = 1.5;                         end
% function val =   q0;    val = 0.8;                         end

% function val = kmin;    val = 1.1;                         end
% function val = kmax;    val = 1.6;                         end
% function val =   q0;    val = 1.05;                         end

function val = zero;    val = 0.0;                         end
function val = one;     val = 1.0;                         end
function val = two;     val = 2.0;                         end
function val = three;   val = 3.0;                         end
function val = four;    val = 4.0;                         end
function val = five;    val = 5.0;                         end
function val = half;    val = one/two;                     end

function val = gamma;   val = 1.4;                         end
function val = gm1;     val = gamma - one;                 end
function val = gp1;     val = gamma + one;                 end
function val = xg;      val = one / gamma;                 end
function val = xgm1;    val = one / (gamma - one);         end
function val = xgp1;    val = one / (gamma + one);         end
function val = xg2m1;   val = one / (gamma*gamma - one);   end
function val = gxgm1;   val = gamma / (gamma - one);       end
function val = gxg2m1;  val = gamma / (gamma*gamma - one); end
function val = gm1xgp1; val = gm1/gp1;                     end
function val = gp1xgm1; val = gp1/gm1;                     end