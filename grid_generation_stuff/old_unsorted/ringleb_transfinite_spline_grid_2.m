function [x,y,psi,V] = ringleb_transfinite_spline_grid_2( N_psi, N_V, ...
                                                          psi_min, ...
                                                          psi_max, V_min)

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
    [psi(:,1:N_V2),V(:,1:N_V2)] = ndgrid( linspace(-1           ,1,N_psi), ...
                                          linspace(-1+1/(N_V2-1),1,N_V2));
    for i = 1:N_V2
        x(:,     i) =  F.x(psi(:,N_V2+1-i),V(:,N_V2+1-i));
        y(:,     i) = -F.y(psi(:,N_V2+1-i),V(:,N_V2+1-i));
        x(:,N_V2+i) =  F.x(psi(:,i),V(:,i));
        y(:,N_V2+i) =  F.y(psi(:,i),V(:,i));
    end
end

[psi,V] = grid_to_psiV(x,y,psi_min,V_min);

end






function [psi,V] = grid_to_psiV(X,Y,psi_min,V_min)
[M,N] = size(X);
V = zeros(M,N);
% psi = zeros(M,N);

maxiter = 100;
tol     = 1e-12;
Vguess = 0.5*(V_min + 1/psi_min);
for j = 1:N
    for i = 1:M
        V(i,j) = find_solution_V(X(i,j),Y(i,j),Vguess,maxiter,tol);
    end
end
V(Y<0) = -V(Y<0);
psi = psi_fun(X,V);
end

function b = b_value(V)
b = sqrt(1 - 0.2*V.^2);
end

function rho = density(V)
b = b_value(V);
rho = b.^5;
end

function psi = psi_fun(xp,V)
    b   = b_value(V);
    rho = b.^5;
    L   = 1./b + 1./(3*b.^3) + 1./(5*b.^5) ...
        - 0.5*log((1+b)./(1-b));
    psi = sqrt(0.5./(V.^2) - rho.*(xp-0.5*L));
end

function L = L_value(V)
b = b_value(V);
L = 1./b + 1./(3*b.^3) + 1./(5*b.^5) - 0.5*log((1+b)./(1-b)); 
end

function x = x_fun(psi,V)
rho = density(V);
L   = L_value(V);
x = 1./rho .* (0.5./(V.^2) - psi.^2) + 0.5*L; 
end

function y = y_fun(psi,V)
rho = density(V);
rad = 1 - (V.^2).*(psi.^2);
rad = (rad>=0).*rad + (rad<0)*0;
y = psi./(rho.*V) .* sqrt(rad);
end


function d = dLdb_fun(b)
den = b.^6.*(b.^2-1);
den = den.*(den~=0) + eps(1).*(den==0);
d = 1./(den);
end

function d = dbdV_fun(V)
d = -(V/5)./sqrt(1 - V.^2/5);
end

function d = drhodV_fun(V)
d = -V.*(1 - V.^2/5).*sqrt(1 - V.^2/5);
end

function d = dxdV_fun(psi,V)
rho    = density(V);
b      = b_value(V);
drhodV = drhodV_fun(V);
dLdV   = dLdb_fun(b).*dbdV_fun(V);
d = (-1./rho.^2).*drhodV.*(1./(2*V.^2)) ...
   + (1./rho).*(-1./V.^3) ...
   + (-1./rho.^2).*drhodV.*(-psi.^2) ...
   + (1/2)*dLdV;
end

% function d = dydV_fun(psi,V)
% rho    = density(V);
% drhodV = drhodV_fun(V);
% small = 1e-12;
% d = ( (-1./rho.^2).*drhodV.*(psi./V) + (-psi./(rho.*V.^2)) ) ...
%     .*sqrt(1-V.^2.*psi.^2) ...
%     + (psi./(rho.*V)) ...
%     .*(-V.*psi.^2./sqrt(1-V.^2.*psi.^2 + small));
% end

function d = dydV_fun(psi,V)
rho    = density(V);
drhodV = drhodV_fun(V);
small = 1e-8;
rad = 1-V.^2.*psi.^2;
d = ( (-1./rho.^2).*drhodV.*(psi./V) + (-psi./(rho.*V.^2)) ) ...
    .*sqrt(1-V.^2.*psi.^2) ...
    + (psi./(rho.*V)) ...
    .*(-V.*psi.^2./sqrt(rad));
d(rad<small) = ( y_fun(psi,V+small) - y_fun(psi,V) )/small;
end

function V2 = reparameterize_streamline(psi,V)
N = length(V);
darcL = @(V) sqrt(dxdV_fun(psi,V).^2 + dydV_fun(psi,V).^2);
% arc length approximation
L = integral( @(V) darcL(V), V(1), V(N),'ArrayValued',true );

% arc length parameterization
L1 = V(N) - V(1);
V2 = V;
for j = 2:N-1
    % f2 = @(t) V(j) - V(j-1) - (L1/L)*integral( @(z) darcL(z), V2(j-1), t,...
    %                                            'ArrayValued',true );
    f2 = @(t) V(j) - V(j-1) - (L1/L)*integral( @(z) darcL(z), V2(j-1), t,...
                                               'ArrayValued',true );
    a1 = V(1); % needs to be a larger value than the previous one
    b1 = V(N);  %
    V2(j) = fzero( f2, [a1,b1]);
end

end


function V2 = reparameterize_streamline_spline(psi,V)
N = length(V);

% darcL = @(V) sqrt( dxdV_fun(psi,V).^2 + dydV_fun(psi,V).^2);
X = x_fun(psi,V);
Y = y_fun(psi,V);
% sx = my_cubic_spline(V,X,1.0e99,1.0e99);
% sy = my_cubic_spline(V,Y,1.0e99,1.0e99);
sx = my_linear_spline(V,X);
sy = my_linear_spline(V,Y);

% darcL = @(V) darcL_fun_spline_approx(sx,sy,V);
% arc length approximation
% L = integral( @(V) darcL(V), V(1), V(N),'ArrayValued',true );
L = int_darcL_fun_linear_spline_approx(sx,sy,V(1),V(N));
% arc length parameterization
L1 = V(N) - V(1);
V2 = V;
for j = 2:N-1
    % f2 = @(t) V(j) - V(j-1) - (L1/L)*integral( @(z) darcL(z), V2(j-1), t,...
    %                                            'ArrayValued',true );
    % f2 = @(t) V(j) - V(j-1) - (L1/L)*integral( @(z) darcL(z), V2(j-1), t,...
    %                                            'ArrayValued',true );
    f2 = @(t) V(j) - V(j-1) - (L1/L)*int_darcL_fun_linear_spline_approx(sx,sy,V2(j-1),t);
    a1 = V(1); % needs to be a larger value than the previous one
    b1 = V(N);  %
    V2(j) = fzero( f2, [a1,b1]);
end

end

function d = darcL_fun_spline_approx(sx,sy,V)
d = sqrt( sx.deval(V).^2 + sy.deval(V).^2 );
end

function d = int_darcL_fun_linear_spline_approx(sx,sy,V1,VN)
dx = sx.L_approx(V1,VN);
dy = sy.L_approx(V1,VN);
d = sum( sqrt( dx.^2 + dy.^2 ) );
if V1 > VN
    d = -d;
end
end

% my_cubic_spline(t,XY{i}(:,1),sx1(i),sx2(i));


function F = construct_transfinite_spline_map_full( N_psi, N_V, psi_min, ...
                                                    psi_max, V_min )
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
% V = reparameterize_streamline_spline(psi,V);
V = reparameterize_streamline(psi,V);
XY2(:,1) = x_fun(psi,V);
XY2(:,2) = y_fun(psi,V);
XY2(1,2) = 0;

% inflow / top_boundary:    V = V_min, k = [psi_min,psi_max]
XY3 = zeros( N_psi, 2 );
psi = linspace( psi_max, psi_min, N_psi );
V   = V_min;
XY3(:,1) = x_fun(psi,V);
XY3(:,2) = y_fun(psi,V);

% outer wall/left boundary: V = [V_min,1/psi_max], psi = psi_max
XY4 = zeros( N_V, 2 );
psi = psi_max;
V = linspace( 1/psi_max, V_min, N_V );
V = reparameterize_streamline(psi,V);
% V = reparameterize_streamline_spline(psi,V);
XY4(:,1) = x_fun(psi,V);
XY4(:,2) = y_fun(psi,V);
XY4(1,2) = 0;
L = 1.0e99;
sx1 = [1,0,L,0];
sy1 = [0,1,L,1];

sx2 = [1,L,L,L];
sy2 = [0,L,L,L];
F = transfinite_spline_interpolant_2D_w_slope_local( XY1, XY2, XY3, XY4, ...
                                                     sx1, sx2, sy1, sy2 );

end

function V = find_solution_V(xp,yp,V0,maxiter,tol)
% initial guess
% Vm1 = 0.5*(V_min + 1/psi_min);
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
    den = 2*rho*sqrt((xp-0.5*L)^2 + yp^2);
    V = sqrt(1/den);
end

function F = transfinite_spline_interpolant_2D_w_slope_local( ...
                     XY1, XY2, XY3, XY4, sx1, sx2, sy1, sy2 )
G = struct();
XY = {XY1,XY2,XY3,XY4};
for i = 1:4
    N = size(XY{i},1);
    t = linspace(-1,1,N);
    G(i).G = struct();
    G(i).G.sx = my_cubic_spline(t,XY{i}(:,1),sx1(i),sx2(i));
    G(i).G.sy = my_cubic_spline(t,XY{i}(:,2),sy1(i),sy2(i));
    G(i).G.x = @(t) G(i).G.sx.eval(t);
    G(i).G.y = @(t) G(i).G.sy.eval(t);
end

F = transfinite_interpolant(G(1).G,G(2).G,G(3).G,G(4).G);

end