function zeta = moretti_map_atan2(z,mu,n)
zm  = z-mu;
zms = z+conj(mu);
d = abs( (zm)./(zms) ).*exp( 1i*(1./n)*( atan2(imag(zm),real(zm)) - atan2(imag(zms),real(zms)) ) );
zeta = (1+d)./(1-d);
end