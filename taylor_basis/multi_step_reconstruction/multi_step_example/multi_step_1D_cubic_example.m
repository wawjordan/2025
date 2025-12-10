%% Coupled Least Squares in 1D (v2: started 05/22/2024)
% from "A high order interpolation for the finite volume method: 
%       The coupled least squares reconstruction" (2018)
% https://doi.org/10.1016/j.compfluid.2018.09.009
clc; clear; close all;

%% Define the mesh
% 1D periodic grid for now
rng('default')
N  = 33;
Nc = N - 1;
x0 = 0; xn = 1;
x  = rand_grd(x0,xn,N,0.4);
xc = ( x(2:end) + x(1:end-1) )/2;
S = calc_grid_moments(x,xc);

%% exact function



% need a periodic function for comparison
f1 = @(x) sin(2*pi*x) + 0.5*cos(10*pi*x);
% f1 = @(x) sinc(10*x)-sin(x)+log(x+2);
% f1 = @(x) (x-0.1).*(x-0.5).*(x-0.9);

%% exact averages
u_avg = get_averages(f1,x);
u_e   = construct_piecewise_const(x,u_avg);

%% compute reconstruction
coefs = get_reconstruct_coefs(S,u_avg);
u_p     = construct_piecewise_poly(coefs,x,xc);

% integral average of the polynomial (should match the exact averages given)
u_p_avg = construct_piecewise_const(x,get_averages(u_p,x));

% error
E1 = @(x) u_p(x) - f1(x);


%% plot
fig2 = figure(2);
fplot(E1, [x0+0.01,xn-0.01],'k')

fig1 = figure(1);
hold on
fplot(f1, [x0+0.01,xn-0.01],'k-','LineWidth',1)
legend_string = {'$f(x)$'};
legend(legend_string,'Interpreter','latex','location','northwest','Box','off')

fplot(u_e,[x0,xn],'k')
legend_string = [legend_string,'$\frac{1}{\vert x_{i+1} - x_i \vert }\int_{x_i}^{x_{i+1}}f(x)dx$'];
legend(legend_string,'Interpreter','latex','Box','off')

fplot(u_p, [x0,xn],'r','LineWidth',1)
legend_string = [legend_string,'$p(x)$'];
legend(legend_string,'Interpreter','latex','Box','off')

fplot(u_p_avg,[x0,xn],'r--')
legend_string = [legend_string,'$\frac{1}{\vert x_{i+1} - x_i \vert }\int_{x_i}^{x_{i+1}}p(x)dx$'];
legend(legend_string,'Interpreter','latex','Box','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coefs = get_reconstruct_coefs(S,u_avg)

sigma1    = step1(u_avg,S.cab1);
sigma_sum = step2(S.cab2,sigma1);
psi       = step3(S.cab2,S.cab3,sigma1);
theta     = step4(S.a2,sigma_sum,psi);
sigma2    = step5(S.a1,S.b1,theta,psi,sigma1);

Nc = length(u_avg);
coefs = zeros(Nc,4);
coefs(:,1) = (1/6)*psi;
coefs(:,2) = (1/2)*theta;
coefs(:,3) = sigma2;
coefs(:,4) = u_avg - (1/24)*theta.*S.tau.^2;

end

function sigma1 = step1(u,cab1)

Nc = length(u);

sigma1 = zeros(Nc,1);

for alpha = 1:Nc
    idx1 = get_stencil(alpha,Nc);
    for beta = 1:3
        idx2 = get_stencil(idx1(beta),Nc);
        sigma1(alpha) = sigma1(alpha) + cab1(alpha,beta)*( u(idx2(2)) - u(alpha) );
    end
end

end

function sigma_sum = step2(cab2,sigma1)

Nc = length(sigma1);

sigma_sum = zeros(Nc,1);

for alpha = 1:Nc
    idx1 = get_stencil(alpha,Nc);
    for beta = 1:3
        idx2 = get_stencil(idx1(beta),Nc);
        sigma_sum(alpha) = sigma_sum(alpha) + cab2(alpha,beta)*( sigma1(idx2(2)) - sigma1(alpha) );
    end
end

end


function psi = step3(cab2,cab3,sigma1)

Nc = length(sigma1);

psi = zeros(Nc,1);

for alpha = 1:Nc
    idx1 = get_stencil(alpha,Nc);
    for beta = 1:3
        idx2 = get_stencil(idx1(beta),Nc);

        tmp_sum1 = 0;
        for gamma = 1:3
            idx3 = get_stencil(idx2(gamma),Nc);
            tmp_sum1 = tmp_sum1 + cab2(idx2(2),gamma)*( sigma1(idx3(2)) - sigma1(idx2(2)) );
        end

        tmp_sum2 = 0;
        for delta = 1:3
            idx3 = get_stencil(idx2(delta),Nc);
            tmp_sum2 = tmp_sum2 + cab2(alpha,delta)*( sigma1(idx3(2)) - sigma1(alpha) );
        end
        psi(alpha) = psi(alpha) + cab3(alpha,beta)*(tmp_sum1-tmp_sum2);
    end
end

end

function theta = step4(a2,sigma_sum,psi)
Nc = length(sigma_sum);
theta = zeros(Nc,1);
for alpha = 1:Nc
    theta(alpha) = -a2(alpha)*psi(alpha) + sigma_sum(alpha);
end
end


function sigma2 = step5(a1,b1,theta,psi,sigma1)
Nc = length(theta);
sigma2 = zeros(Nc,1);
for alpha = 1:Nc
    sigma2(alpha) = - a1(alpha)*theta(alpha) ...
                    - b1(alpha)*psi(alpha)   ...
                    + sigma1(alpha);
end
end


function S = calc_grid_moments(x,xc)
% x and xc are cell nodes and centers
% N is the number of moments for each stencil
% (size of the stencil is 3)

N = 4;

Nc = length(xc);

% pre-allocate
tau  = zeros(Nc,1);
h    = zeros(N,Nc,3);
H    = zeros(N,Nc);
cab1 = zeros(Nc,3);
cab2 = zeros(Nc,3);
cab3 = zeros(Nc,3);
a1   = zeros(Nc,1);
b1   = zeros(Nc,1);
a2   = zeros(Nc,1);
z    = zeros(N,Nc,3);


%%
for alpha = 1:Nc
    tau(alpha) = x( alpha + 1 ) - x( alpha );
end

%%
for beta = 1:3
    for alpha = 1:Nc
        idx1 = get_stencil(alpha,Nc);
        idx2 = get_stencil(idx1(beta),Nc);
        h(1,alpha,beta) = xc( idx2(2) ) - xc( alpha );
    end
end

alpha = 1;
beta  = 1;
idx1 = get_stencil(alpha,Nc);
idx2 = get_stencil(idx1(1),Nc);
h(1,alpha,beta) = xc( idx2(2) ) - xc( alpha ) - ( x(idx2(2)+1) - x(alpha) );

alpha = Nc;
beta  = 3;
idx1 = get_stencil(alpha,Nc);
idx2 = get_stencil(idx1(3),Nc);
h(1,alpha,beta) = xc( idx2(2) ) - xc( alpha ) - ( x(idx2(2)) - x(alpha+1) );

for beta = 1:3
    for alpha = 1:Nc
        for k = 2:N
            h(k,alpha,beta) = h(1,alpha,beta)^k;
        end
    end
end

%%
for alpha = 1:Nc
    for k = 1:N
        H(k,alpha) = sum(h(k,alpha,:));
    end
end

%%
for beta = 1:3
    for alpha = 1:Nc
        cab1(alpha,beta) = h(1,alpha,beta)/sum(h(1,alpha,:).^2);
    end
end

%%
for alpha = 1:Nc
    a1(alpha) = (1/2)*H(3,alpha)/H(2,alpha);
    b1(alpha) = (1/6)*H(4,alpha)/H(2,alpha);
    idx1 = get_stencil(alpha,Nc);
    for beta = 1:3
        idx2 = get_stencil(idx1(beta),Nc);
        a1(alpha) = a1(alpha) + (1/24)*cab1(alpha,beta)*( tau(idx2(2))^2 - tau(alpha)^2 );
        b1(alpha) = b1(alpha) + (1/24)*cab1(alpha,beta)*tau(idx2(2))^2*h(1,alpha,beta);
    end
end

%%
for alpha = 1:Nc
    a2(alpha) = a1(alpha);
    den       = 1;
    idx1 = get_stencil(alpha,Nc);
    for beta = 1:3
        idx2 = get_stencil(idx1(beta),Nc);
        a2(alpha) = a2(alpha) + cab1(alpha,beta)*( b1(idx2(2)) - b1(alpha) + a1(idx2(2))*h(1,alpha,beta) );
        den       = den + cab1(alpha,beta)*( a1(idx2(2)) - a1(alpha) );
    end
    a2(alpha) = a2(alpha)/den;
    cab2(alpha,:) = cab1(alpha,:)/den;
end

%%
for alpha = 1:Nc
    den       = 1;
    idx1 = get_stencil(alpha,Nc);
    for delta = 1:3
        idx2 = get_stencil(idx1(beta),Nc);
        den       = den + cab1(alpha,delta)*( a2(idx2(2)) - a2(alpha) );
    end
    cab3(alpha,:) = cab1(alpha,:)/den;
end


%%
for alpha = 1:Nc
    idx = get_stencil(alpha,Nc);
    % for each cell in the stencil
    for beta = 1:3
        x1 = x( idx(beta)     ) - xc( alpha ); % lower bound
        x2 = x( idx(beta) + 1 ) - xc( alpha ); % upper bound
        % integrate the moment
        for k = 1:N
            z(k,beta,alpha) = (1/(k+1)/tau(beta))*( x2^k - x1^k );
        end
    end
end

S = struct();
S.tau  = tau;
S.h    = h;
S.H    = H;
S.cab1 = cab1;
S.cab2 = cab2;
S.cab3 = cab3;
S.a1   = a1;
S.b1   = b1;
S.a2   = a2;
S.z    = z;

end

function wi = get_stencil(i,n)
% compact stencil
wi = mod([i-1,i,i+1]-1,n) + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = rand_grd(x0,xn,n,C)
% x0 - start value
% xn - stop  value
%  n - number of points
%  C - irregularity constant

% kappa = 1 + C*rand(n-1,1);
kappa = 1 + C*(1-2*rand(n-1,1));
h = (xn - x0)/sum(kappa);

x = zeros(n,1);
x(1) = x0;
x(2:n) = x0 + cumsum(kappa*h);

end

function u_e = construct_piecewise_const(x,u_avg)
    idx_fun = @(xq) discretize(xq,x,'IncludedEdge','left');
    u_e = @(xq) arrayfun(@(xq)u_avg(idx_fun(xq)),xq);
end

function u_p = construct_piecewise_poly(coefs,x,xc)

idx_fun = @(xq) discretize(xq,x,'IncludedEdge','left');
u_eval = @(xq,i) eval_poly( coefs(i,:), xc(i),xq);

u_p = @(xq) arrayfun(@(xq)u_eval(xq,idx_fun(xq)),xq);
end

function u_avg = get_averages(fun,x)
Nc = length(x)-1;
u_avg = zeros(Nc,1);
for i = 1:Nc
    vol = abs(x(i+1) - x(i));
    u_avg(i) = (1/vol)*integral(fun,x(i),x(i+1));
end
end

function u_eval = eval_poly(coefs,xc,x)
u_eval = 0;
xtmp = x - xc;
Ncoef = length(coefs);
for k = 1:Ncoef-1
    u_eval = u_eval + coefs(k)*xtmp^(Ncoef-k);
end
u_eval = u_eval + coefs(Ncoef);
end