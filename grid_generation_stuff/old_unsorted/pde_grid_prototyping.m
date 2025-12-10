%% 1D grid spacing differential methods
clc; clear; close all;

% define the intervals




phi{1} = @(x) sin(pi*x);
phi{2} = @(x) x.*sin(pi*x).^2;
phi{3} = @(x) tanh(x);
phi{4} = @(x) 1.1 + cos(2*pi*x);
phi{5} = @(x) exp(-20*(x-0.5).^2);

labels = cellfun(@func2str,phi,'UniformOutput',false);

M = length(phi);

xmin = -1; xmax = 1;
N = 101;


x = zeros(N,M);

for j = 1:M
    x(:,j) = weight_grid(xmin,xmax,N,phi{j});
end

figure(1)
hold on
for j = 1:M
    plot(x(:,j));
    % plot(diff( x(:,j) ));
end
legend(labels)

figure(2)
hold on
for j = 1:M
    plot(phi{j}(linspace(0,1,N)));
end
legend(labels)


% function x = weight_grid(xmin,xmax,N,phi)
% xi = linspace(0,1,N);
% xi_c = 0.5*( xi(1:N-1) + xi(2:N) );
% abc = 1./phi(xi_c);
% di = zeros(N-2,1);
% di(1) = abc(1)*xmin;
% di(N-2) = abc(N-1)*xmax;
% x = zeros(N,1);
% x(1) = xmin;
% x(N) = xmax;
% x(2:N-1) = local_tridiag(-abc(1:N-2),abc(1:N-2) + abc(2:N-1),-abc(2:N-1),di);
% end



function x = local_tridiag(a,b,c,d)
N = length(b);
x = zeros(N,1);
for i = 2:N
    w = a(i)/b(i-1);
    b(i) = b(i)-w*c(i-1);
    d(i) = d(i)-w*d(i-1);
end
x(N) = d(N)/b(N);
for i = N-1:-1:1
    x(i) = (d(i)-c(i)*x(i+1))/b(i);
end
% M = zeros(N);
% M(N+1:N+1:N*N-1) = a(2:N);
% M(1:N+1:N*N) = b;
% M(2:N+1:N*(N-1)) = c(1:N-1);
end


% ai = zeros(m-2,1); bi = zeros(m-2,1); ci = zeros(m-2,1);
% 
% for i = 2:m-1
%     phi_m12 = phi( 0.5*( xi(i-1) + xi(i) ) );
%     ai(i-1) = 1/phi_m12;
% end
% for i = 2:m-1
%     phi_m12 = phi( 0.5*( xi(i-1) + xi(i) ) );
%     phi_p12 = phi( 0.5*( xi(i+1) + xi(i) ) );
%     bi(i-1) = 1/phi_m12 + 1/phi_p12;
% end
% for i = 2:m-1
%     phi_p12 = phi( 0.5*( xi(i+1) + xi(i) ) );
%     ci(i-1) = 1/phi_p12;
% end
% di = zeros(m-2,1);
% di(1) = ai(1)*a;
% di(m-2) = ci(m-2)*b;
% 
% x = zeros(m,1);
% x(1) = a;
% x(m) = b;
% x(2:m-1) = local_tridiag(-ai,bi,-ci,di);
% 
% xi_c = 0.5*( xi(1:m-1) + xi(2:m) );
% phi_c = phi(xi_c);
% aci = 1./phi_c;
% bi = aci(2:m-1) + aci(1:m-2);