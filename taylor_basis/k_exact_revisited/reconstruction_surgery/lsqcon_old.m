function x = lsqcon_old(A, c, B, d)
%lsqcon_old    Constrained least squares.
%   x = lsqcon_old(A, b, B, d) solves the constrained least squares problem
%   min|Ax-b|_2 st Bx = d.

[Q, U, R, S] = gqr_old(A, B, 0);
[~,n] = size(A);  p = size(S,2);
i1 = 1:p;         i2 = (p+1):n;
y1 = S(i1,i1).'\d;
y2 = R(i2,i2).'\(U(i2,:)*c-R(i1,i2).'*y1);
x = Q*[y1 ; y2];

end