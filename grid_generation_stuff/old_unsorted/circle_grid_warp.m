%% Circle grid warping 04/10/24
clc; clear; close all;

% l          = 1;
% epsilon    = 0.1;
% kappa      = -0.0;
% tau        = deg2rad(0.174532925199428);
l          = 1;
epsilon    = 0.3;
kappa      = -0.0;
tau        = pi/4;

scale      = false;
a0         = l*sqrt( (1+epsilon)^2 + kappa^2 );
zeta_c     = -l*epsilon + l*kappa*1i;
beta       = asin(kappa/a0);


% N_theta = 65;
% N_r = 49;
N_theta = 4097;
N_r = 3073;

r_outer = 100;
dTE = 0.5;

bump_fun = @(x,b) exp(-1./(b.^2-x.^2)).*(abs(x)<b)*exp(1/b^2);

mag1 = 2;

b = 2;
wt_fun = @(x) (mag1-1)*bump_fun(x*2*b-b,b) + 1;
% wt_fun = @(x) bump_fun(x*2-1,1)+0.1;
theta_space = @(xmin,xmax,N) weight_grid(xmin,xmax,N,wt_fun);
% theta_space = @(xmin,xmax,N) build_and_eval_tanh_fun(xmin,xmax,N,dTE*(xmax-xmin)/N,dTE*(xmax-xmin)/N,@linspace);
% theta_space = @(xmin,xmax,N) linspace(xmin,xmax,N);



c = 4;

% theta_fun = @(theta) (mag-1)*( 1 + cos(theta) )/2 + 1;
% theta_fun = @(theta) (mag-1)*( 1 + cos(theta) )/2 + 1;
% theta_fun = @(theta) (mag-1)*( 2-exp(-0.5*(theta-pi).^2) )/2 + 1;
% theta_fun = @(theta) (mag-1)*bump_fun(theta-pi,tau/2) + 1;



% theta_fun = @(theta) (mag-1)*(1-bump_fun(theta-pi,pi-tau/10000)) + 1;
% expansion_ratio = @(theta) (c.*theta_fun(theta).*r_outer).^(1/N_r);
% expansion_ratio = @(theta) (c.*r_outer).^(1/N_r);

mag2 = 2;

base_expansion_ratio = (c*r_outer).^(1/N_r);
te_expansion_ratio = (c*mag2*r_outer).^(1/N_r);
expansion_ratio = @(theta) ( ( te_expansion_ratio - base_expansion_ratio).*(1-bump_fun(theta-pi,pi)) + 1 )*base_expansion_ratio;



theta1 = theta_space(0,2*pi,N_theta);
theta2 = linspace(0,2*pi,N_theta);
R2 = linspace( a0, 4*r_outer, N_r );
R     = zeros(N_theta,N_r);
THETA = zeros(N_theta,N_r);
for i = 1:N_theta
    R(i,:) = mygeomspace( a0, 4*r_outer, expansion_ratio(theta1(i)), N_r );
    
    % tmp = R(i,:)./R2;
    % tmp = (tmp-min(tmp))/(max(tmp)-min(tmp));
    % THETA(i,:) = (1-tmp).*theta2(i) + tmp.*theta1(i);
    THETA(i,:) = theta1(i);
end

[Z,ZETA] = kt_grid_function(R,THETA,l,epsilon,kappa,tau,scale);

skip = 4;
Z1 = Z(1:2^skip:end,1:2^skip:end);
ZETA1 = ZETA(1:2^skip:end,1:2^skip:end);


xi = real(ZETA1);
eta = imag(ZETA1);

x = real(Z1);
y = imag(Z1);

% figure(1)
% hold on;
% plot(xi,eta,'k');
% plot(xi.',eta.','k');
% axis equal

figure(2)
hold on;
plot(x,y,'k');
plot(x.',y.','k');
xlim([-2,2])
ylim([-2,2])
axis equal



function t = build_and_eval_tanh_fun(x1,x2,N,d1,d2,spacing)
% F = tanh_stretching_function(x1,x2,d1,d2,N,spacing);

F = tanh_stretching_function(x1,x2,d1,d2,N);
t = F(spacing(x1,x2,N));

end