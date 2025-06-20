function x = lsqcon_alt_old(A, b, B, d)
% compute the QR factorization B^T=QR
[Q,R] = qr(B.');
% m1 = size(A,1);
n1 = size(A,2);
m2 = size(B,1);

% solve R(1:m2,1:m2)^T*y = d for y
y = R(1:m2,1:m2).' \ d;

% A = A*Q
A = A*Q;

% find z so || A(:,m2+1:n1)*z - ( b - A(:,1:m2)*y ) ||_2 is minimized
z = A(:,m2+1:n1) \ ( b - A(:,1:m2)*y );

x = Q(:,1:m2)*y + Q(:,m2+1:n1)*z;

end