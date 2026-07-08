function x = bilinear_map(A,B,C,D,uv)
% X(u,v) = A + (B-A)⋅u + (D-A)⋅v + (A-B+C-D)⋅u⋅v
% P = A + (B-A).*u;
% Q = C + (C-D).*u;
% x = P + (Q-P).*v;
x = A + (B-A)*uv(1) + (D-A)*uv(2) + (A-B+C-D)*uv(1)*uv(2);
end