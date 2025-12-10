function [X2,Y2,R] = winslow_eqns_fd(X,Y,max_iter,tol)
X2 = X; Y2 = Y;
[Ni,Nj] = size(X);
R = zeros(max_iter,2);

[X2,Y2,RESX,RESY] = update_xy(X2,Y2,Ni,Nj);
resinit(1) = sqrt(sum(RESX.^2,'all')/(Ni*Nj));
resinit(2) = sqrt(sum(RESY.^2,'all')/(Ni*Nj));
R(1,:) = resinit;
if all( resinit < 4*eps(1) )
    return
end
for i = 2:max_iter
    [X2,Y2,RESX,RESY] = update_xy(X2,Y2,Ni,Nj);
    resnorm(1) = sqrt(sum(RESX.^2,'all')/(Ni*Nj));
    resnorm(2) = sqrt(sum(RESY.^2,'all')/(Ni*Nj));
    R(i,:) = resnorm./resinit;
    fprintf('%d  %e  %e\n',i,R(i,1),R(i,2))
    if all(R(i,:) < tol )
        fprintf('Converged!\n');
        return
    end
end

end


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

function [RESX,RESY] = residual_calc(X,Y,Ni,Nj)
RESX = zeros(Ni-2,Nj-2);
RESY = zeros(Ni-2,Nj-2);
% constant j
for j = 2:Nj-1
    for i = 2:Ni-1
        g11 = g11_diff(X(i-1:i+1,j),Y(i-1:i+1,j));
        g22 = g22_diff(X(i,j-1:j+1),Y(i,j-1:j+1));
        g12 = g12_diff(X(i-1:i+1,j-1:j+1),Y(i-1:i+1,j-1:j+1));
        AC = g22;
        B = ( 2*g22 + 2*g11 );
        D = g11*( X(i  ,j+1) + X(i  ,j-1) ) ...
             - 2*g12*( X(i+1,j+1) + X(i-1,j-1) ...
                     - X(i-1,j+1) - X(i+1,j-1) )/4;
        E = g11*( Y(i  ,j+1) + Y(i  ,j-1) ) ...
             - 2*g12*( Y(i+1,j+1) + Y(i-1,j-1) ...
                     - Y(i-1,j+1) - Y(i+1,j-1) )/4;
        RESX(i-1,j-1) = -AC*X(i-1,j) + B*X(i,j) - AC*X(i+1,j) - D;
        RESY(i-1,j-1) = -AC*Y(i-1,j) + B*Y(i,j) - AC*Y(i+1,j) - E;
    end
end

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

function [X,Y,RESX,RESY] = update_xy(X,Y,Ni,Nj)
% i - implicit sweep
% for j = 2:Nj-1
%     [A,B,C,D,E] = approx_coefs_j(Ni,X(:,j-1:j+1),Y(:,j-1:j+1));
%     X(2:Ni-1,j) = local_tridiag(A,B,C,D);
%     Y(2:Ni-1,j) = local_tridiag(A,B,C,E);
% end

for j = 2:Nj-1
    [A,B,C,D,E] = approx_coefs_j(Ni,X(:,j-1:j+1),Y(:,j-1:j+1));
    X(2:Ni-1,j) = local_tridiag(A,B,C,D);
    Y(2:Ni-1,j) = local_tridiag(A,B,C,E);
end


% D = zeros(Nj-2,1);

% j - implicit sweep
for i = 2:Ni-1
    [A,B,C,D,E] = approx_coefs_i(Nj,X(i-1:i+1,:),Y(i-1:i+1,:));
    X(i,2:Nj-1) = local_tridiag(A,B,C,D);
    Y(i,2:Nj-1) = local_tridiag(A,B,C,E);
end
[RESX,RESY] = residual_calc(X,Y,Ni,Nj);
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