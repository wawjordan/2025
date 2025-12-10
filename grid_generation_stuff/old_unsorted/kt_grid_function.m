function [Z,ZETA] = kt_grid_function(R,THETA,l,epsilon,kappa,tau,scale)

% radius of circle in zeta plane
a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
n = 2 - tau/pi;
zeta_c = -l*epsilon + l*kappa*1i;
beta = asin(kappa/a0);

ZETA = R.*exp(1i*(THETA - beta)) + zeta_c;
Z = KT_zeta_to_z(ZETA,n,l);

if scale
    [c,zLE,~] = get_airfoil_chord(l,epsilon,kappa,tau);
    Z = (Z - zLE)/c;
end

end

function [c,zLE,zTE] = get_airfoil_chord(l,epsilon,kappa,tau)
% radius of circle in zeta plane
a0 = l*sqrt( (1+epsilon)^2 + kappa^2 );
beta = asin(kappa/a0);
n = 2 - tau/pi;
zeta_c = -l*epsilon + l*kappa*1i;
z_fun = @(s,theta)s*real(KT_zeta_to_z(a0.*exp(1i*(theta - beta)) + zeta_c,n,l));
theta_LE = fminbnd(@(theta)z_fun( 1,theta),0,2*pi);
theta_TE = fminbnd(@(theta)z_fun(-1,theta),0,2*pi);

zeta_LE = a0.*exp(1i*(theta_LE - beta)) + zeta_c;
zeta_TE = a0.*exp(1i*(theta_TE - beta)) + zeta_c;
zLE = KT_zeta_to_z(zeta_LE,n,l);
zTE = KT_zeta_to_z(zeta_TE,n,l);
c = ( abs(zLE) + abs(zTE) );
end

function z = KT_zeta_to_z(zeta,n,l)
zeta_p = (zeta+l).^n;
zeta_m = (zeta-l).^n;

z = n*l*( (zeta_p + zeta_m)./(zeta_p - zeta_m) );

end

% function zeta = KT_z_to_zeta(z,n,l)
% z_p = (z+n*l);
% z_m = (z-n*l);
% 
% zeta = -l*( (z_m./z_p).^(1/n) + 1) ./ ( (z_m./z_p).^(1/n) - 1);
% 
% end

% function d = diff_KT_zeta_to_z(zeta,n,l)
% zeta_frac = ( (zeta - l)./(zeta + l) ).^n;
% 
% d = 4*(n*l)^2*zeta_frac./( (zeta.^2-1).*(1-zeta_frac).^2 );
% end
