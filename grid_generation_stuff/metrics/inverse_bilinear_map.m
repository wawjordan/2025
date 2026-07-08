function [uv,stat] = inverse_bilinear_map(A,B,C,D,x)
% https://iquilezles.org/articles/ibilinear/
tol = 1e-12;
E = B(:)-A(:);
F = D(:)-A(:);
G = A(:)-B(:)+C(:)-D(:);
H = x(:)-A(:);
k2 = G(1)*F(2) - G(2)*F(1);
k1 = E(1)*F(2) - E(2)*F(1) + H(1)*G(2) - H(2)*G(1);
k0 = H(1)*E(2) - H(2)*E(1);

% k2 = cross2d(G,F);
% k1 = cross2d(E,F) + cross2d(H,G);
% k0 = cross2d(H,E);
stat = 0;
uv = [0;0];
if (abs(k2)<tol)
    uv(1) = ( H(1)*k1 + F(1)*k0)/(E(1)*k1-G(1)*k0);
    uv(2) = -k0/k1;
else
    discr = k1*k1 - 4*k0*k2;
    if (discr<0) % outside of [0,1]
        stat = 1;
        return
    end
    discr = sqrt(discr);
    xk2 = 0.5/k2;
    uv(2) = (-k1 - discr)*xk2;
    uv(1) = ( H(1) - F(1)*uv(2) )/( E(1) + G(1)*uv(2) );

    if ( uv(1)<0 || uv(1)>1 || uv(2)<0 || uv(2)>1 ) % choose the other root
        uv(2) = (-k1 + discr)*xk2;
        uv(1) = ( H(1) - F(1)*uv(2) ) / ( E(1) + G(1)*uv(2) );
    end
end

end
% 
% function val = cross2d(v1,v2)
% tmp = cross([v1(:);0],[v2(:);0]);
% val = tmp(3);
% end