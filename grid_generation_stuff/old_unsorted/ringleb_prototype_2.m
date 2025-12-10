%% Ringleb Flow prototyping - from "I do like CFD"
clc; clear; close all;

% spacing1 = {@(xmin,xmax,N)mygeomspace( xmin, xmax, 1.00, N ),@(xmin,xmax,N)mygeomspace( xmin, xmax, 0.99, N )};
N_theta = 129;
N_psi   = 33;

% psi_min = 0.7;
% psi_max = 1.5;
% V_min   = 0.4;

% subsonic
psi_min = 2;
psi_max = 6;
V_min   = 0.12;

% supersonic
% psi_min = 0.6;
% psi_max = 0.8;
% V_min   = 1.2;

% [X,Y,psi,V] = basic_grid_1(N_theta,N_psi,psi_min,psi_max,V_min);
% [X,Y,psi,V] = transfinite_spline_grid(N_psi,N_theta,psi_min,psi_max,V_min);
% [X2,Y2,psi2,V2] = basic_grid_2(N_theta,N_psi,psi_min,psi_max,V_min);

% [X1,Y1,~,~] = basic_grid_1_full(N_theta,N_psi,psi_min,psi_max,V_min);
% [X2,Y2,~,~] = basic_grid_2_full(N_theta,N_psi,psi_min,psi_max,V_min);
[X,Y,psi,V] = transfinite_spline_grid_full(N_psi,N_theta,psi_min,psi_max,V_min);

% hold on; axis equal;
% plot(X,Y,'k'); plot(X.',Y.','k');
% plot(X1,Y1,'r'); plot(X1.',Y1.','r');
% plot(X2,Y2,'g'); plot(X2.',Y2.','g');


% [Xc,Yc,psic,Vc] = center_grid_approx(X,Y,V);
S = get_soln(X,V);
M = abs(V)./sqrt(S(:,:,5)./S(:,:,1));
% Sc = get_soln(Xc,Vc);
% Mc = abs(Vc)./sqrt(Sc(:,:,5)./Sc(:,:,1));
hold on; axis equal;
contourf(X,Y,M,200,'EdgeColor','none');
colorbar
plot(X,Y,'k'); plot(X.',Y.','k');



function [X,Y,psi,V] = basic_grid_1(N_theta,N_psi,psi_min,psi_max,V_min)
psi_0 = linspace( psi_min, psi_max, N_psi );
V   = zeros(N_theta,N_psi);
psi = zeros(N_theta,N_psi);

for i = 1:N_theta
    psi(i,:) = psi_0;
end

for j = 1:N_psi
    theta = linspace(asin(psi_0(j)*V_min),pi/2,N_theta);
    V(:,j) = sin(theta)/psi_0(j);
end

V(N_theta,:) = one./psi_0;

X = x_fun(psi,V);
Y = y_fun(psi,V);

end

function [X,Y,psi,V] = basic_grid_1_full(N_theta,N_psi,psi_min,psi_max,V_min)

psi_0 = linspace( psi_min, psi_max, N_psi );
V   = zeros(N_theta,N_psi);
psi = zeros(N_theta,N_psi);

for i = 1:N_theta
    psi(i,:) = psi_0;
end

for j = 1:N_psi
    theta_min = asin(psi_0(j)*V_min);
    theta_max = pi - asin(psi_0(j)*V_min);
    theta = linspace(theta_min,theta_max,N_theta);
    V(:,j) = sin(theta)/psi_0(j);
end

if mod(N_theta,2) == 1
    V((N_theta+1)/2,:) = one./psi_0;
    V((N_theta+1)/2+1:N_theta,:) = -V((N_theta+1)/2+1:N_theta,:);
else
    V(N_theta/2+1:N_theta,:) = -V(N_theta/2+1:N_theta,:);
end

X = x_fun(psi,V);
Y = y_fun(psi,V);

end

function [X,Y,psi,V] = basic_grid_2(N_theta,N_psi,psi_min,psi_max,V_min)
psi_0 = linspace( psi_min, psi_max, N_psi );
V   = zeros(N_theta,N_psi);
psi = zeros(N_theta,N_psi);

for i = 1:N_theta
    psi(i,:) = psi_0;
end

for j = 1:N_psi
    theta_min = asin(V_min*psi_0(j));
    phi_start = 1./(psi_0(j)*V_min)  *cos(theta_min);
    phi_end   = 0;
    phi = linspace(phi_start,phi_end,N_theta);
    for i = 1:N_theta
        k = 1/psi_0(j);
        V(i,j) = sqrt(k^2/(1+phi(i)^2));
    end
end

X = x_fun(psi,V);
Y = y_fun(psi,V);

end

function [X,Y,psi,V] = basic_grid_2_full(N_theta,N_psi,psi_min,psi_max,V_min)
psi_0 = linspace( psi_min, psi_max, N_psi );
V   = zeros(N_theta,N_psi);
psi = zeros(N_theta,N_psi);

for i = 1:N_theta
    psi(i,:) = psi_0;
end

for j = 1:N_psi
    theta_min = asin(V_min*psi_0(j));
    phi_start = 1./(psi_0(j)*V_min)  *cos(theta_min);
    phi_end   = -phi_start;
    phi = linspace(phi_start,phi_end,N_theta);
    for i = 1:N_theta
        k = 1/psi_0(j);
        V(i,j) = sqrt(k^2/(1+phi(i)^2));
    end
end

if mod(N_theta,2) == 1
    V((N_theta+1)/2,:) = one./psi_0;
    V((N_theta+1)/2+1:N_theta,:) = -V((N_theta+1)/2+1:N_theta,:);
else
    V(N_theta/2+1:N_theta,:) = -V(N_theta/2+1:N_theta,:);
end

X = x_fun(psi,V);
Y = y_fun(psi,V);

end


function [Xc,Yc,psic,Vc] = center_grid_approx(X,Y,V)
[M,N] = size(V);
Xc = zeros(M-1,N-1);
Yc = zeros(M-1,N-1);
Vc = V(1:M-1,1:N-1);
% psic = zeros(N-1,M-1);

maxiter = 100;
tol     = 1e-12;

for j = 1:N-1
    for i = 1:M-1
        Xc(i,j) = ( X(i,j) + X(i+1,j) + X(i,j+1) + X(i+1,j+1) )/4;
        Yc(i,j) = ( Y(i,j) + Y(i+1,j) + Y(i,j+1) + Y(i+1,j+1) )/4;
        Vc(i,j) = find_solution_V(Xc(i,j),Yc(i,j),Vc(i,j),maxiter,tol);
    end
end
Vc(Yc<0) = -Vc(Yc<0);
% Vm1 = half*(V_min + 1/psi_min);
% Vc = arrayfun(@(xp,yp,V0)find_solution_V(xp,yp,V0,maxiter,tol),Xc,Yc,Vc);
psic = psi_fun(Xc,Vc);
end


function [psi,V] = grid_to_psiV(X,Y,psi_min,V_min)
[M,N] = size(X);
V = zeros(M,N);
% psi = zeros(M,N);

maxiter = 100;
tol     = 1e-12;
Vguess = half*(V_min + 1/psi_min);
for j = 1:N
    for i = 1:M
        V(i,j) = find_solution_V(X(i,j),Y(i,j),Vguess,maxiter,tol);
    end
end
V(Y<0) = -V(Y<0);
psi = psi_fun(X,V);
end





function b = b_value(V)
b = sqrt(one - 0.2*V.^2);
end

function rho = density(V)
b = b_value(V);
rho = b.^5;
end

function psi = psi_fun(xp,V)
    b   = b_value(V);
    rho = b.^5;
    L   = one./b + one./(three*b.^3) + one./(five*b.^5) ...
        - half*log((one+b)./(one-b));
    psi = sqrt(half./(V.^2) - rho.*(xp-half*L));
end

function L = L_value(V)
b = b_value(V);
L = one./b + one./(three*b.^3) + one./(five*b.^5) - half*log((one+b)./(one-b)); 
end

function x = x_fun(psi,V)
rho = density(V);
L   = L_value(V);
x = one./rho .* (half./(V.^2) - psi.^2) + half*L; 
end

function y = y_fun(psi,V)
rho = density(V);
rad = one - (V.^2).*(psi.^2);
rad = (rad>=0).*rad + (rad<0)*0;
y = psi./(rho.*V) .* sqrt(rad);
end


function d = dLdb_fun(b)
den = b.^6.*(b.^2-one);
den = den.*(den~=0) + eps(one).*(den==0);
d = one./(den);
end

function d = drhodb_fun(b)
d = 5*b.^4;
end

function d = dbdV_fun(V)
d = -(V/5)./sqrt(one - V.^2/5);
end

function d = drhodV_fun(V)
d = -V.*(one - V.^2/5).*sqrt(one - V.^2/5);
% b = b_value(V);
% d2 = -V.*b.^3;
% d3 = drhodb_fun(b).*dbdV_fun(V);
end

function d = dxdV_fun(psi,V)
rho    = density(V);
b      = b_value(V);
drhodV = drhodV_fun(V);
dLdV   = dLdb_fun(b).*dbdV_fun(V);
% d      = ( (psi.^2./rho.^2 - one./(rho.*V).^2).*drhodV - one./(rho.*V.^3) + dLdV)/2;
d = (-1./rho.^2).*drhodV.*(1./(2*V.^2)) ...
   + (1./rho).*(-1./V.^3) ...
   + (-1./rho.^2).*drhodV.*(-psi.^2) ...
   + (1/2)*dLdV;
% plot(V,d2,'g')
end

function d = dydV_fun(psi,V)
rho    = density(V);
drhodV = drhodV_fun(V);
small = 1e-6;
d = ( (-1./rho.^2).*drhodV.*(psi./V) + (-psi./(rho.*V.^2)) ).*sqrt(1-V.^2.*psi.^2) ...
    + (psi./(rho.*V)).*(-V.*psi.^2./sqrt(1-V.^2.*psi.^2+ small));
end

function d = dxdpsi_fun(psi,V)
rho = density(V);
d = - two.*psi./rho;
end

function d = dydpsi_fun(psi,V)
rho = density(V);
rad = one - (V.^2).*(psi.^2);
rad = (rad>=0).*rad + (rad<0)*0;
sqrt(rad);
d = sqrt(rad)./(rho.*V) - (psi.^4)./(rho.*V)./sqrt(rad); 
end

function V2 = reparameterize_streamline(psi,V)
N = length(V);

darcL = @(V) sqrt(dxdV_fun(psi,V).^2 + dydV_fun(psi,V).^2);
% arc length approximation
L = integral( @(V) darcL(V), V(1), V(N),'ArrayValued',true,'RelTol',1e-14 );

% arc length parameterization
L1 = V(N) - V(1);
V2 = V;
for j = 2:N-1
    f2 = @(t) V(j) - V(j-1) - (L1/L)*integral( @(z) darcL(z), V2(j-1), t,'ArrayValued',true );
    a1 = V(1); % needs to be a larger value than the previous one
    b1 = V(N);  %
    % x1 = linspace( a1, b1, 10 );
    % hold on
    % plot(x1,arrayfun(f2,x1))
    V2(j) = fzero( f2, [a1,b1]);
end

end

function V2 = reparameterize_streamline_alt(psi,V)
N = length(V);
step = 1.0e-8;
darcL = @(V) sqrt( ( ( x_fun(psi,V+step) - x_fun(psi,V) )/step ).^2 ...
                 + ( ( y_fun(psi,V+step) - y_fun(psi,V) )/step ).^2 );
% arc length approximation
L = integral( @(V) darcL(V), V(1), V(N),'ArrayValued',true );

% arc length parameterization
L1 = V(N) - V(1);
V2 = V;
for j = 2:N-1
    f2 = @(t) V(j) - V(j-1) - (L1/L)*integral( @(z) darcL(z), V2(j-1), t,'ArrayValued',true );
    a1 = V(1); % needs to be a larger value than the previous one
    b1 = V(N);  %
    % x1 = linspace( a1, b1, 10 );
    % hold on
    % plot(x1,arrayfun(f2,x1))
    V2(j) = fzero( f2, [a1,b1]);
end

end

function psi2 = reparameterize_potential_line(psi,V)
N = length(psi);

darcL = @(psi) sqrt(dxdpsi_fun(psi,V).^2 + dydpsi_fun(psi,V).^2);
% arc length approximation
L = integral( @(psi) darcL(psi), psi(1), psi(N),'ArrayValued',true );

% arc length parameterization
L1 = psi(N) - psi(1);
psi2 = psi;
for j = 2:N-1
    f2 = @(t) psi(j) - psi(j-1) - (L1/L)*integral( @(z) darcL(z), psi2(j-1), t,'ArrayValued',true );
    a1 = psi(1);
    b1 = psi(N);  %
    % x1 = linspace( a1, b1, 10 );
    % hold on
    % plot(x1,arrayfun(f2,x1))
    psi2(j) = fzero( f2, [a1,b1]);
end

end

function F = construct_transfinite_spline_map(N_psi,N_V,psi_min,psi_max,V_min)
% symmetry/bottom boundary: V = 1/psi, psi = [psi_min,psi_max]
XY1 = zeros( N_psi, 2 );
psi = linspace( psi_max, psi_min, N_psi );
V   = 1./psi;
XY1(:,1) = linspace( x_fun(psi(1),V(1)), x_fun(psi(N_psi),V(N_psi)), N_psi );
XY1(:,2) = 0;

% inner wall/right boundary: - V = [V_min,1/psi_min], psi = psi_min
XY2 = zeros( N_V, 2 );
psi = psi_min;
V = linspace( 1/psi_min, V_min, N_V );
V = reparameterize_streamline(psi,V);
XY2(:,1) = x_fun(psi,V);
XY2(:,2) = y_fun(psi,V);
XY2(1,2) = 0;

% inflow / top_boundary:    V = V_min, k = [psi_min,psi_max]
XY3 = zeros( N_psi, 2 );
psi = linspace( psi_max, psi_min, N_psi );
V   = V_min;
% psi = reparameterize_potential_line(psi,V);
XY3(:,1) = x_fun(psi,V);
XY3(:,2) = y_fun(psi,V);

% outer wall/left boundary: V = [V_min,1/psi_max], psi = psi_max
XY4 = zeros( N_V, 2 );
psi = psi_max;
% V = linspace( V_min, 1/psi_max, N_psi );
V = linspace( 1/psi_max, V_min, N_V );
V = reparameterize_streamline(psi,V);
XY4(:,1) = x_fun(psi,V);
XY4(:,2) = y_fun(psi,V);
XY4(1,2) = 0;
F = transfinite_spline_interpolant_2D(XY1,XY2,XY3,XY4);
end

% function F = construct_transfinite_spline_map_full(N_psi,N_V,psi_min,psi_max,V_min)
% if mod(N_V,2) == 1
%     N_V2 = (N_V+1)/2;
%     XY2 = zeros( N_V, 2 );
%     XY4 = zeros( N_V, 2 );
% else
%     N_V2 = N_V/2 + 1;
%     XY2 = zeros( N_V+1, 2 );
%     XY4 = zeros( N_V+1, 2 );
% end
% % inner wall/right boundary: - V = [V_min,1/psi_min], psi = psi_min
% 
% psi = psi_min;
% V = linspace( V_min, 1/psi_min, N_V2 );
% V = reparameterize_streamline(psi,V);
% 
% for i = 1:N_V2-1
%     XY2(i,1) =  x_fun( psi, V(i) );
%     XY2(i,2) = -y_fun( psi, V(i) );
% end
% 
% i = N_V2;
% 
% XY2(i,1) = x_fun(psi,V(i));
% XY2(i,2) = 0;
% 
% for i = 1:N_V2-1
%     XY2(N_V2+i,1) = x_fun( psi, V(N_V2-i) );
%     XY2(N_V2+i,2) = y_fun( psi, V(N_V2-i) );
% end
% 
% % outer wall/left boundary: - V = [V_min,1/psi_max], psi = psi_max
% psi = psi_max;
% V = linspace( V_min, 1/psi_max, N_V2 );
% V = reparameterize_streamline(psi,V);
% 
% 
% for i = 1:N_V2-1
%     XY4(i,1) =  x_fun( psi, V(i) );
%     XY4(i,2) = -y_fun( psi, V(i) );
% end
% 
% i = N_V2;
% 
% XY4(i,1) = x_fun(psi,V(i));
% XY4(i,2) = 0;
% 
% for i = 1:N_V2-1
%     XY4(N_V2+i,1) = x_fun( psi, V(N_V2-i) );
%     XY4(N_V2+i,2) = y_fun( psi, V(N_V2-i) );
% end
% 
% % outflow / bottom_boundary:    V = V_min, k = [psi_min,psi_max]
% XY1 = zeros( N_psi, 2 );
% psi = linspace( psi_max, psi_min, N_psi );
% V   = V_min;
% % psi = reparameterize_potential_line(psi,V);
% XY1(:,1) = x_fun(psi,V);
% XY1(:,2) = -y_fun(psi,V);
% 
% % inflow / top_boundary:    V = V_min, k = [psi_min,psi_max]
% XY3 = zeros( N_psi, 2 );
% psi = linspace( psi_max, psi_min, N_psi );
% V   = V_min;
% % psi = reparameterize_potential_line(psi,V);
% XY3(:,1) = x_fun(psi,V);
% XY3(:,2) = y_fun(psi,V);
% 
% % hold on;
% % plot(XY1(:,1),XY1(:,2),'r')
% % plot(XY1(1,1),XY1(1,2),'r.')
% % plot(XY2(:,1),XY2(:,2),'b')
% % plot(XY2(1,1),XY2(1,2),'b.')
% % plot(XY3(:,1),XY3(:,2),'g')
% % plot(XY3(1,1),XY3(1,2),'g.')
% % plot(XY4(:,1),XY4(:,2),'c')
% % plot(XY4(1,1),XY4(1,2),'c.')
% 
% F = transfinite_spline_interpolant_2D(XY1,XY2,XY3,XY4);
% 
% end

function F = construct_transfinite_spline_map_full(N_psi,N_V,psi_min,psi_max,V_min)
% symmetry/bottom boundary: V = 1/psi, psi = [psi_min,psi_max]
XY1 = zeros( N_psi, 2 );
psi = linspace( psi_max, psi_min, N_psi );
V   = 1./psi;
XY1(:,1) = linspace( x_fun(psi(1),V(1)), x_fun(psi(N_psi),V(N_psi)), N_psi );
XY1(:,2) = 0;

% inner wall/right boundary: - V = [V_min,1/psi_min], psi = psi_min
XY2 = zeros( N_V, 2 );
psi = psi_min;

V0 = linspace( 1/psi_min, V_min, N_V );
V = reparameterize_streamline(psi,V0);
% V1 = reparameterize_streamline_alt(psi,V0);

XY2(:,1) = x_fun(psi,V);
XY2(:,2) = y_fun(psi,V);
XY2(1,2) = 0;

% inflow / top_boundary:    V = V_min, k = [psi_min,psi_max]
XY3 = zeros( N_psi, 2 );
psi = linspace( psi_max, psi_min, N_psi );
V   = V_min;
% psi = reparameterize_potential_line(psi,V);
XY3(:,1) = x_fun(psi,V);
XY3(:,2) = y_fun(psi,V);

% outer wall/left boundary: V = [V_min,1/psi_max], psi = psi_max
XY4 = zeros( N_V, 2 );
psi = psi_max;
% V = linspace( V_min, 1/psi_max, N_psi );
V = linspace( 1/psi_max, V_min, N_V );
V = reparameterize_streamline(psi,V);
XY4(:,1) = x_fun(psi,V);
XY4(:,2) = y_fun(psi,V);
XY4(1,2) = 0;
F = transfinite_spline_interpolant_2D(XY1,XY2,XY3,XY4);
L = 1.0e99;
sx1 = [1,0,L,0];
sy1 = [0,1,L,1];

sx2 = [1,L,L,L];
sy2 = [0,L,L,L];
F = transfinite_spline_interpolant_2D_w_slope(XY1,XY2,XY3,XY4,sx1,sx2,sy1,sy2);

end


function [x,y,psi,V] = transfinite_spline_grid(N_psi,N_V,psi_min,psi_max,V_min)
F = construct_transfinite_spline_map(N_psi,N_V,psi_min,psi_max,V_min);
[psi,V] = ndgrid(linspace(-1,1,N_psi),linspace(-1,1,N_V));
x = F.x(psi,V);
y = F.y(psi,V);

[psi,V] = grid_to_psiV(x,y,psi_min,V_min);
end

function [x,y,psi,V] = transfinite_spline_grid_full(N_psi,N_V,psi_min,psi_max,V_min)
F = construct_transfinite_spline_map_full(N_psi,N_V,psi_min,psi_max,V_min);
x = zeros(N_psi,N_V);
y = zeros(N_psi,N_V);
if mod(N_V,2) == 1
    N_V2 = (N_V+1)/2;
    [psi,V] = ndgrid(linspace(-1,1,N_psi),linspace(-1,1,N_V2));
    for i = 1:N_V2
        x(:,     i) =  F.x(psi(:,N_V2+1-i),V(:,N_V2+1-i));
        y(:,     i) = -F.y(psi(:,N_V2+1-i),V(:,N_V2+1-i));
        x(:,N_V2-1+i) =  F.x(psi(:,i),V(:,i));
        y(:,N_V2-1+i) =  F.y(psi(:,i),V(:,i));
    end
else
    N_V2 = N_V/2;
    [psi(:,1:N_V2),V(:,1:N_V2)] = ndgrid(linspace(-1,1,N_psi),linspace(-1+1/(N_V2-1),1,N_V2));
    for i = 1:N_V2
        x(:,     i) =  F.x(psi(:,N_V2+1-i),V(:,N_V2+1-i));
        y(:,     i) = -F.y(psi(:,N_V2+1-i),V(:,N_V2+1-i));
        x(:,N_V2+i) =  F.x(psi(:,i),V(:,i));
        y(:,N_V2+i) =  F.y(psi(:,i),V(:,i));
    end
end

[psi,V] = grid_to_psiV(x,y,psi_min,V_min);

end

function state = exact_soln_nondim(xp,V)
    b   = b_value(V);
    rho = b.^5;
    p   = b.^7;
    L   = one./b + one./(three*b.^3) + one./(five*b.^5) ...
        - half*log((one+b)./(one-b));
    psi   = sqrt(half./V.^2 - rho.*(xp-half*L));
    arg = psi*V;
    arg = arg*(abs(arg)<=1) + 1*(arg>1) - 1*(arg<-1);
    theta = asin(arg);
    u     = -V*cos(theta);
    v     = -V*sin(theta);
    state = zero;
    state(1) = rho;
    state(2) = u;
    state(3) = v;
    state(5) = p;
end


function S = get_soln(X,V)
[M,N] = size(X);

S = zeros(M,N,5);
for j = 1:N
    for i = 1:M
        S(i,j,:) = exact_soln_nondim(X(i,j),V(i,j));
    end
end
        
end




function V = find_solution_V(xp,yp,V0,maxiter,tol)
% initial guess
% Vm1 = half*(V_min + 1/psi_min);
V = Vguess(xp,yp,V0);
err0 = abs(V-V0);
V0 = V;
for i = 1:maxiter
    V = Vguess(xp,yp,V0);
    R = abs(V-V0)/err0;
    if R < tol
        break
    end
    V0 = V;
end

if (R > tol)
    warning('tol not met in find_solution_V')
end

end

function V = Vguess(xp,yp,V)
    rho = density(V);
    L   = L_value(V);
    den = two*rho*sqrt((xp-half*L)^2 + yp^2);
    V = sqrt(one/den);
end

function val = zero;    val = 0.0;                         end
function val = one;     val = 1.0;                         end
function val = two;     val = 2.0;                         end
function val = three;   val = 3.0;                         end
function val = four;    val = 4.0;                         end
function val = five;    val = 5.0;                         end
function val = half;    val = one/two;                     end

function val = gamma;   val = 1.4;                         end