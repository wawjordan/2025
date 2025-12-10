%% 2D Discretized Winslow Equations (from Basic Structured Grid Gen.) 11/01/23
clc; clear; close all;

kmin = 0.7;
kmax = 1.5;
q0   = 0.4;
spacing = { @(xmin,xmax,N)mygeomspace( xmin, xmax, 1.0, N ),...
             @(xmin,xmax,N)mygeomspace( xmin, xmax, 1.01, N ) };
% spacing = { @linspace, @linspace };
Nq = 129;
Nk = 65;

% [X,Y,~,~] = ringleb_streamline_grid(kmin,kmax,q0,Nq,Nk,spacing);
[X,Y,~,~] = ringleb_transfinite_spline_grid(kmin,kmax,q0,Nq,Nk,spacing);
% [X,Y,~,~] = ringleb_transfinite_grid(kmin,kmax,q0,Nq,Nk,spacing);
X1 = X; Y1 = Y;
for i = 1:100
    [X1,Y1] = update_xy(X1,Y1,Nq,Nk);
end

hold on; axis equal
plot(X,Y,'k'); plot(X.',Y.','k');
plot(X1,Y1,'r--'); plot(X1.',Y1.','r--');

function g11 = g11_diff(x,y)
i = 2; j = 1;
g11 = ( ( x(i+1,j) - x(i-1,j) )/2 )^2 + ( ( y(i+1,j) - y(i-1,j) )/2 )^2;
end

function g22 = g22_diff(x,y)
i = 1; j = 2;
g22 = ( ( x(i,j+1) - x(i,j-1) )/2 )^2 + ( ( y(i,j+1) - y(i,j-1) )/2 )^2;
end

function g12 = g12_diff(x,y)
i = 2; j = 2;
g12 = ( ( x(i+1,j) - x(i-1,j) )/2 )*( ( x(i,j+1) - x(i,j-1) )/2 ) ...
    + ( ( y(i+1,j) - y(i-1,j) )/2 )*( ( y(i,j+1) - y(i,j-1) )/2 );
end

function [A,B,C,D,E] = approx_coefs_i(Nj,X,Y)
% constant j
AC = zeros(Nj,1);
B = zeros(Nj-2,1);
D = zeros(Nj-2,1);
E = zeros(Nj-2,1);
i = 2;
for j = 2:Nj-1
    g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
    g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
    g12 = g12_diff(X(i-1:i+1,j-1:j+1),Y(i-1:i+1,j-1:j+1));

    AC(j) = g11;
    B(j-1) = ( 2*g22 + 2*g11 );
    D(j-1) = g22*( X(i+1,j  ) + X(i-1,j  ) ) ...
         - 2*g12*( X(i+1,j+1) + X(i-1,j-1) ...
                 - X(i-1,j+1) - X(i+1,j-1) )/4;
    E(j-1) = g22*( Y(i+1,j  ) + Y(i-1,j  ) ) ...
         - 2*g12*( Y(i+1,j+1) + Y(i-1,j-1) ...
                 - Y(i-1,j+1) - Y(i+1,j-1) )/4;
end

for j = [1,Nj]
    g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
    AC(j) = g11;
end


A = -AC(2:Nj-1);
C = -AC(2:Nj-1);
D(1)    = D(1)    + AC(2)   *X(2, 1);
D(Nj-2) = D(Nj-2) + AC(Nj-1)*X(2,Nj);
% 
E(1)    = E(1)    + AC(2)   *Y(2, 1);
E(Nj-2) = E(Nj-2) + AC(Nj-1)*Y(2,Nj);

end


function [A,B,C,D,E] = approx_coefs_j(Ni,X,Y)
% constant j
AC = zeros(Ni,1);
B = zeros(Ni-2,1);
D = zeros(Ni-2,1);
E = zeros(Ni-2,1);
j = 2;
for i = 2:Ni-1
    g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
    g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
    g12 = g12_diff(X(i-1:i+1,j-1:j+1),Y(i-1:i+1,j-1:j+1));

    AC(i) = g22;
    B(i-1) = ( 2*g22 + 2*g11 );
    D(i-1) = g11*( X(i  ,j+1) + X(i  ,j-1) ) ...
         - 2*g12*( X(i+1,j+1) + X(i-1,j-1) ...
                 - X(i-1,j+1) - X(i+1,j-1) )/4;
    E(i-1) = g11*( Y(i  ,j+1) + Y(i  ,j-1) ) ...
         - 2*g12*( Y(i+1,j+1) + Y(i-1,j-1) ...
                 - Y(i-1,j+1) - Y(i+1,j-1) )/4;
end

for i = [1,Ni]
    g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
    AC(i) = g22;
end


A = -AC(2:Ni-1);
C = -AC(2:Ni-1);
D(1)    = D(1)    + AC(2)   *X( 1,2);
D(Ni-2) = D(Ni-2) + AC(Ni-1)*X(Ni,2);
% 
E(1)    = E(1)    + AC(2)   *Y( 1,2);
E(Ni-2) = E(Ni-2) + AC(Ni-1)*Y(Ni,2);
end


% function [A,B,C,D,E] = approx_coefs_i(Nj,X,Y)
% % constant j
% AC = zeros(Nj,1);
% B = zeros(Nj-2,1);
% D = zeros(Nj-2,1);
% E = zeros(Nj-2,1);
% i = 2;
% for j = 2:Nj-1
%     g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
%     g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
%     g12 = g12_diff(X(i-1:i+1,j-1:j+1),Y(i-1:i+1,j-1:j+1));
% 
%     AC(j) = g11;
%     B(j-1) = ( 2*g22 + 2*g11 );
%     D(j-1) = g22*( X(i  ,j+1) + X(i  ,j-1) ) ...
%          - 2*g12*( X(i+1,j+1) + X(i-1,j-1) ...
%                  - X(i-1,j+1) - X(i+1,j-1) )/4;
%     E(j-1) = g22*( Y(i  ,j+1) + Y(i  ,j-1) ) ...
%          - 2*g12*( Y(i+1,j+1) + Y(i-1,j-1) ...
%                  - Y(i-1,j+1) - Y(i+1,j-1) )/4;
% end
% 
% for j = [1,Nj]
%     g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
%     AC(j) = g11;
% end
% 
% 
% A = -AC(2:Nj-1);
% C = -AC(2:Nj-1);
% D(1)    = D(1)    + AC(2)   *X(2, 1);
% D(Nj-2) = D(Nj-2) + AC(Nj-1)*X(2,Nj);
% % 
% E(1)    = E(1)    + AC(2)   *Y(2, 1);
% E(Nj-2) = E(Nj-2) + AC(Nj-1)*Y(2,Nj);
% 
% end
% 
% 
% function [A,B,C,D,E] = approx_coefs_j(Ni,X,Y)
% % constant j
% AC = zeros(Ni,1);
% B = zeros(Ni-2,1);
% D = zeros(Ni-2,1);
% E = zeros(Ni-2,1);
% j = 2;
% for i = 2:Ni-1
%     g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
%     g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
%     g12 = g12_diff(X(i-1:i+1,j-1:j+1),Y(i-1:i+1,j-1:j+1));
% 
%     AC(i) = g22;
%     B(i-1) = ( 2*g22 + 2*g11 );
%     D(i-1) = g11*( X(i+1,j  ) + X(i-1,j  ) ) ...
%          - 2*g12*( X(i+1,j+1) + X(i-1,j-1) ...
%                  - X(i-1,j+1) - X(i+1,j-1) )/4;
%     E(i-1) = g11*( Y(i+1,j  ) + Y(i-1,j  ) ) ...
%          - 2*g12*( Y(i+1,j+1) + Y(i-1,j-1) ...
%                  - Y(i-1,j+1) - Y(i+1,j-1) )/4;
% end
% 
% for i = [1,Ni]
%     g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
%     AC(i) = g22;
% end
% 
% 
% A = -AC(2:Ni-1);
% C = -AC(2:Ni-1);
% D(1)    = D(1)    + AC(2)   *X( 1,2);
% D(Ni-2) = D(Ni-2) + AC(Ni-1)*X(Ni,2);
% % 
% E(1)    = E(1)    + AC(2)   *Y( 1,2);
% E(Ni-2) = E(Ni-2) + AC(Ni-1)*Y(Ni,2);
% end


function [X,Y] = update_xy(X,Y,Ni,Nj)
XC = X;
YC = Y;
% i - implicit sweep
for j = 2:Nj-1
    [A,B,C,D,E] = approx_coefs_j(Ni,X(:,j-1:j+1),Y(:,j-1:j+1));
    X(2:Ni-1,j) = local_tridiag(A,B,C,D);
    Y(2:Ni-1,j) = local_tridiag(A,B,C,E);
end

% j - implicit sweep
for i = 2:Ni-1
    [A,B,C,D,E] = approx_coefs_i(Nj,X(i-1:i+1,:),Y(i-1:i+1,:));
    X(i,2:Nj-1) = local_tridiag(A,B,C,D);
    Y(i,2:Nj-1) = local_tridiag(A,B,C,E);
end

end

function x = local_tridiag(a,b,c,d)
N = length(b);

% M = zeros(N);
% M(N+1:N+1:N*N-1) = a(2:N);
% M(1:N+1:N*N) = b;
% M(2:N+1:N*(N-1)) = c(1:N-1);

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
end