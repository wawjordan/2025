%% KT airfoil grid prototype (03/29/2024)
clc; clear; close all;

inputs = struct();
inputs.l          = 1;
inputs.epsilon    = 0.1;
inputs.kappa      = -0.0;
% inputs.tau        = deg2rad(0.174532925199428);
inputs.tau        = 0.0; %0.174532925199428;


% inputs.l          = 1;
% inputs.epsilon    = 0.1;
% inputs.kappa      = 0.2;
% inputs.tau        = 0.1;

rhoinf = 1.0;
pinf   = 100000.0;
vinf   = 75.0;
alpha  = 0.0;
alpha_rad = deg2rad(alpha);

w_soln = @(x,y) airfoil_velocity(x,y,inputs.l,inputs.epsilon,inputs.kappa,inputs.tau,alpha_rad,vinf);

w_soln2 = @(z) airfoil_velocity_zeta(z,inputs.l,inputs.epsilon,inputs.kappa,inputs.tau,alpha_rad,vinf);

u_soln = @(x,y)  real( w_soln(x,y) );
v_soln = @(x,y) -imag( w_soln(x,y) );

u_soln2 = @(z)  real( w_soln2(z) );
v_soln2 = @(z) -imag( w_soln2(z) );

% soln = @(x,y) airfoil_solution(x,y,rhoinf,vinf,pinf,alpha_rad,gamma,radius,mu,n,l);

% inputs.l          = 1;
% inputs.epsilon    = 0.1;
% inputs.kappa      = 0.2;
% inputs.tau        = 0.1;

% dtheta = 2*pi/(N_theta-1);
% expansion_ratio = (2+dtheta)/(2-dtheta);

% in the limit as N_theta goes to infinity
%  ((2+b)/(2-b))^(alpha*N_theta) = exp(2*pi*alpha)
% where b = 2pi/N_theta
% and alpha = N_r/N_theta

% if we have a desired outer radius, and N_theta, then
% r_outer = 100;
% N_r = log(r_outer)/log(expansion_ratio);

% if we have a desired outer radius, N_theta, and N_r, then
r_outer = 100;
% N_theta = 1025;
% N_r = 769;
N_theta = 4097;
N_r = 3073;

% N_theta = 65;
% N_r = 49;


% get the chord of the transformed airfoil
[c,~,~] = get_airfoil_chord(inputs.l,inputs.epsilon,inputs.kappa,inputs.tau);
[c2,~,~] = get_airfoil_chord(inputs.l/c,inputs.epsilon,inputs.kappa,inputs.tau);

expansion_ratio = (c*r_outer)^(1/N_r);
% expansion_ratio = (c*4*r_outer)^(1/N_r);
inputs.r_outer = r_outer;

% dLE = 1;
dTE = 0.2;
% d_r0 = 0.01;
% d_r1 = 3;
inputs.spacing{1} = @linspace;
% inputs.spacing{1} = @(xmin,xmax,N) build_and_eval_tanh_fun(xmin,xmax,N,d_r0*(xmax-xmin)/N,d_r1*(xmax-xmin)/N,@linspace);
% inputs.spacing{2} = @(xmin,xmax,N) build_and_eval_dist_fun(xmin,xmax,N,dLE*(xmax-xmin)/N,dTE*(xmax-xmin)/N,@linspace);

% inputs.spacing{1} = @(xmin,xmax,N)mygeomspace( xmin, xmax, expansion_ratio, N );
inputs.spacing{2} = @linspace;
% inputs.spacing{2} = @(xmin,xmax,N) build_and_eval_tanh_fun(xmin,xmax,N,dTE*(xmax-xmin)/N,dTE*(xmax-xmin)/N,@linspace);
inputs.scale      = false;

zeta_c = -inputs.l*inputs.epsilon + inputs.l*inputs.kappa*1i;

[GRID,ZETA] = kt_grid(N_r,N_theta,inputs);

Z1 = GRID.x+1i*GRID.y;
Z2 = KT_zeta_to_z( KT_z_to_zeta(Z1,2 - inputs.tau/pi,inputs.l) - zeta_c/c ,2 - inputs.tau/pi,inputs.l);
ERR = abs(Z2-Z1).^2;

% ZETA2 = KT_z_to_zeta(KT_zeta_to_z(ZETA,2 - inputs.tau/pi,inputs.l),2 - inputs.tau/pi,inputs.l);
% ERR = abs(ZETA2-ZETA).^2;
%%
skip = 5;
X1 = GRID.x(1:2^skip:end,1:2^skip:end);
Y1 = GRID.y(1:2^skip:end,1:2^skip:end);
ZETA1 = ZETA(1:2^skip:end,1:2^skip:end);

u = u_soln(X1,Y1);
v = v_soln(X1,Y1);
p = pinf + 0.5*rhoinf*(vinf^2 - u.^2 - v.^2);

load("u_compare.mat")
a_ref = 374.165738677;
u = u/a_ref;

u2 = u_soln2(ZETA1-zeta_c/c);
u2 = u2/a_ref;

%%
figure(1)
hold on;
contourf(X1,Y1,u,201,'EdgeColor','None');
colormap('hsv')
colorbar;
clim([0,0.25])
plot(X1,Y1,'k');
plot(X1.',Y1.','k');
axis equal
xlim([-2.2,1.8])
ylim([-2,2])
% ylim([-0.5,0.5])

% figure(2)
% hold on;
% contourf(X1,Y1,u_compare,201,'EdgeColor','None');
% colormap('hsv')
% colorbar;
% clim([0,0.25])
% plot(X1,Y1,'k');
% plot(X1.',Y1.','k');
% axis equal
% xlim([-2.2,-1.8])
% ylim([-0.5,0.5])


% u2 = u_soln2(ZETA(1:2^skip:end,1:2^skip:end));
% v2 = v_soln2(ZETA(1:2^skip:end,1:2^skip:end));
% figure(2)
% hold on;
% contourf(X1,Y1,u2,21,'EdgeColor','None');
% colormap('jet')
% colorbar;
% % plot(X1,Y1,'k');
% % plot(X1.',Y1.','k');
% axis equal
% xlim([-3,3])
% ylim([-3,3])



function [GRID,ZETA] = kt_grid(N_r,N_t,inputs)
l       = inputs.l;
epsilon = inputs.epsilon;
kappa   = inputs.kappa;
tau     = inputs.tau;
r_outer = inputs.r_outer;

% radius of circle in zeta plane
a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
n = 2 - tau/pi;
zeta_c = -l*epsilon + l*kappa*1i;
beta = asin(kappa/a0);

% zeta_LE = a0.*exp(1i*(pi - beta)) + zeta_c;
% zeta_TE = a0.*exp(1i*(0 - beta)) + zeta_c;
% zLE1 = KT_zeta_to_z(zeta_LE,n,l);
% zTE1 = KT_zeta_to_z(zeta_TE,n,l);

[c,zLE,~] = get_airfoil_chord(l,epsilon,kappa,tau);

r     = inputs.spacing{1}(a0,c*r_outer,N_r);
theta = inputs.spacing{2}(0,2*pi,N_t);
[R,THETA] = ndgrid(r,theta);
ZETA = R.*exp(1i*(THETA - beta)) + zeta_c;
Z = KT_zeta_to_z(ZETA,n,l);

% r1     = inputs.spacing{1}(a0,c*r_outer,N_r);
% theta = inputs.spacing{2}(0,2*pi,N_t);
% ZETA = zeros(N_r,N_t);
% for j = 1:N_t
%     r2 = reparameterize_curve_r(r1,theta(j),n,l,beta,zeta_c);
%     ZETA(:,j) = r2.*exp(1i*(theta(j) - beta)) + zeta_c;
% end
% Z = KT_zeta_to_z(ZETA,n,l);





if inputs.scale
    Z = (Z - zLE)/c;
end


GRID.imax = N_r;
GRID.jmax = N_t;
GRID.x = real(Z);
GRID.y = imag(Z);
end

function [c,zLE,zTE] = get_airfoil_chord(l,epsilon,kappa,tau)
% radius of circle in zeta plane
a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
beta = asin(kappa/a0);
n = 2 - tau/pi;
zeta_c = -l*epsilon + l*kappa*1i;
z_fun = @(s,theta)s*real(KT_zeta_to_z(a0.*exp(1i*(theta - beta)) + zeta_c,n,l));
theta_LE = fminbnd(@(theta)z_fun( 1,theta),0,2*pi);
theta_TE = fminbnd(@(theta)z_fun(-1,theta),0,2*pi);

zeta_LE = a0.*exp(1i*(theta_LE - beta)) + zeta_c;
zeta_TE = a0.*exp(1i*(theta_TE - beta)) + zeta_c;
zLE = KT_zeta_to_z(zeta_LE,n,l);
zTE = KT_zeta_to_z(zeta_TE,n,l);
c = ( abs(zLE) + abs(zTE) );
end

function z = KT_zeta_to_z(zeta,n,l)
zeta_p = (zeta+l).^n;
zeta_m = (zeta-l).^n;

z = n*l*( (zeta_p + zeta_m)./(zeta_p - zeta_m) );

end

function zeta = KT_z_to_zeta(z,n,l)
z_p = (z+n*l);
z_m = (z-n*l);

zeta = -l*( (z_m./z_p).^(1/n) + 1) ./ ( (z_m./z_p).^(1/n) - 1);

end

function d = diff_KT_zeta_to_z1(zeta,n,l)
zeta_plus  = (zeta + l);
zeta_minus = (zeta - l);

d = 4*(n*l)^2*( (zeta_plus.*zeta_minus).^(n-1) )./ ...
               ( zeta_plus.^n - zeta_minus.^n ).^2;
end

function d = diff_KT_zeta_to_z_wikipedia(zeta,n,l)
zeta_plus  = (1 + 1./zeta).^n;
zeta_minus = (1 - 1./zeta).^n;

d = 4*(n*l)^2*( zeta_plus.*zeta_minus )./ ...
          ( zeta_plus.^n - zeta_minus.^n ).^2;
end

function d = diff_KT_zeta_to_z(zeta,n,l)
zeta_frac = ( (zeta - l)./(zeta + l) ).^n;

d = 4*(n*l)^2*zeta_frac./( (zeta.^2-1).*(1-zeta_frac).^2 );
end


function w = cylinder_solution(z, z1, vinf, alpha, gamma, radius)

alpha_im = 1i*alpha;
gamma_im = 1i*gamma;
w = vinf*exp(-alpha_im) + gamma_im./( 2*pi*(z-z1) ) ...                       &
  - vinf*radius^2*exp(alpha_im)./(z-z1).^2;
end

% function q = airfoil_solution(x,y,rhoinf,vinf,pinf,alpha,gamma,radius,mu,n,l)
% 
% q = zeros(5,1);
% 
% q(1) = rhoinf;
% 
% z = x + 1i*y;
% wtilde = cylinder_solution(z, mu, vinf, alpha, gamma, radius);
% 
% zeta = KT_zeta_to_z(z,n,l);
% dz = diff_KT_zeta_to_z(zeta,n,l);
% 
% w = wtilde./dz;
% 
% q(2) = real(w);
% q(3) = -imag(w);
% % q(5) = pinf + 0.5*rhoinf*(vinf^2 - q(2).^2 - q(3).^2)
% q(5) = pinf + 0.5*rhoinf*(vinf^2 - abs(w).^2);
% 
% end

function w = airfoil_velocity(x,y,l,epsilon,kappa,tau,alpha,vinf)

z = x + 1i*y;

n    = 2 - tau/pi;
zeta = KT_z_to_zeta(z,n,l);
w = airfoil_velocity_zeta(zeta,l,epsilon,kappa,tau,alpha,vinf);
end

function w = airfoil_velocity_zeta(zeta,l,epsilon,kappa,tau,alpha,vinf)
a    = l*sqrt( (1+epsilon)^2 + kappa^2 );
beta = asin(l*kappa/a);
mu   = l*(-epsilon + kappa*1i);
n    = 2 - tau/pi;
dz = diff_KT_zeta_to_z(zeta,n,l);
wtilde = exp(-1i*alpha) + 2*1i*a*sin(alpha+beta)./(zeta-mu) - a^2*exp(1i*alpha)./(zeta-mu).^2;
w = vinf*wtilde./dz;
end

% function [x,y] = KT_airfoil_surface(t,l,tau,epsilon,kappa)
% a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
% n = 2 - tau/pi;
% zeta_c = -l*epsilon + l*kappa*1i;
% beta = asin(kappa/a0);
% zeta = a0.*exp(1i*(t - beta)) + zeta_c;
% z = KT_zeta_to_z(zeta,n,l);
% x = real(z);
% y = imag(z);
% end
% 
% function t = KT_airfoil_surface_inv(x,y,l,tau,epsilon,kappa)
% a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
% n = 2 - tau/pi;
% zeta_c = -l*epsilon + l*kappa*1i;
% beta = asin(kappa/a0);
% 
% z = x + 1i*y;
% zeta = KT_z_to_zeta(z,n,l);
% % zeta = zeta - zeta_c;
% % t = atan2(real(zeta),imag(zeta));
% t = atan(real(zeta)./-imag(zeta));
% 
% t = t + beta;
% t = mod(t-min(t),2*pi);
% end

% function t = build_and_eval_dist_fun(x1,x2,N,dLE,dTE,spacing)
% % F = tanh_stretching_function(x1,x2,d1,d2,N);
% 
% 
% y = [x1,(x1+x2)/2,x2];
% d =  [dTE,dLE,dTE];
% N2 = [floor(N/2),ceil(N/2)];
% 
% 
% F = grid_distribution_function(y,d,N2);
% t = F(spacing(x1,x2,N));
% 
% end
% 
function t = build_and_eval_tanh_fun(x1,x2,N,d1,d2,spacing)
% F = tanh_stretching_function(x1,x2,d1,d2,N,spacing);

F = tanh_stretching_function(x1,x2,d1,d2,N);
t = F(spacing(x1,x2,N));

end

function r2 = reparameterize_curve_r(r1,theta,n,l,beta,zeta_c)
N = length(r1);

circ  = @(r,theta) r.*exp(1i*(theta - beta)) + zeta_c;
darcL = @(r) abs( diff_KT_zeta_to_z(circ(r+1e-12,theta),n,l) );
% arc length approximation
L = integral( @(r) darcL(r), r1(1), r1(N),'ArrayValued',true,'RelTol',1e-14 );

% arc length parameterization
L1 = r1(N) - r1(1);
r2 = r1;
for j = 2:N-1
    f2 = @(t) r2(j) - r1(j-1) - (L1/L)*integral( @(z) darcL(z), r2(j-1), t,'ArrayValued',true );
    a1 = r1(1); % needs to be a larger value than the previous one
    b1 = r2(N);  %
    % x1 = linspace( a1, b1, 10 );
    % hold on
    % plot(x1,arrayfun(f2,x1))
    r2(j) = fzero( f2, [a1,b1]);
end

end
