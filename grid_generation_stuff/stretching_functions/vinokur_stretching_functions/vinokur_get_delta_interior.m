function [ delta, xi1, b ] = vinokur_get_delta_interior(ti,si,refine)
delta = 2*vinokur_sinh_function( 2*si*sqrt(ti*(1-ti)) );
if (refine)
    delta = fzero(@(z) z - 1./sqrt((((cosh(z)-1+1/ti)./sinh(z)).^2-1).*(si*ti)^2),delta);
end
xi1 = log( (1+ti*( exp(delta)-1) )./(1-ti*(1-exp(-delta)) ) )./(2*delta);
b = sinh(xi1*delta);
end

% function delta = refine(ti,si,delta0)
% f = @(delta) delta - eqn(ti,si,delta);
% delta = fzero(f,delta0);
% end
% 
% 
% function delta = eqn(ti,si,delta0)
% delta = 1./sqrt( ( ( ( cosh(delta0) - 1 + 1/ti )./sinh(delta0) ).^2 - 1 ).*(si*ti)^2 );
% end