function [ xi, dxi, d2xi ] = vinokur_interior(t,ti,si,refine)
[ delta, xi1, b ] = vinokur_get_delta_interior(ti,si,refine);
xi = xi1 + asinh( (t/ti - 1)*b )/delta;
dxi = (b/(ti*delta))./sqrt(b^2*(t/ti - 1).^2 + 1);
d2xi = -(b^3*(t/ti - 1))./(ti^2*(b^2*(t/ti - 1).^2 + 1).^(3/2))/delta;
end

% function [ xi, dxi ] = vinokur_interior(t,ti,si)
% delta0 = 2*vinokur_sinh_function( 2*si*sqrt(ti*(1-ti)) );
% delta = refine(ti,si,delta0);
% xi1 = log( (1+ti*( exp(delta)-1) )./(1-ti*(1-exp(-delta)) ) )./(2*delta);
% % xi = xi1 + asinh( (t/ti - 1)*sinh(xi1*delta) )/delta;
% b = sinh(xi1*delta);
% xi = xi1 + asinh( (t/ti - 1)*b )/delta;
% dxi = (b/(ti*delta))./sqrt(b^2*(t/ti - 1).^2 + 1);
% end
% 
% function delta = refine(ti,si,delta0)
% f = @(delta) delta - eqn(ti,si,delta);
% delta = fzero(f,delta0);
% end
% 
% 
% function delta = eqn(ti,si,delta0)
% delta = 1./sqrt( ( ( ( cosh(delta0) - 1 + 1/ti )./sinh(delta0) ).^2 - 1 ).*(si*ti)^2 );
% end