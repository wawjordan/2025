function x = lsqcon_alt2_old(A, b, B, d)
% compute the QR factorization APb=Q[R,f;0,g]
[Q,R,P] = qr(A);

% n1 = size(A,2);
f = Q \ b;

% solve RP^Ty = f  for y
y = R*P.'\f;

% solve (P)(R^T)(K^T) = B^T  for K^T
KT = P*R.' \ B.';

% compute the minimum norm solution of Ku = d -  By
u = KT.' \ (d - B*y);

% solve (R)(P^T)z = u  for z
z = R*P.' \ u;

% set x = y + z
x = y + z;

end