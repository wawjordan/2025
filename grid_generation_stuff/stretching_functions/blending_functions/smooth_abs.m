function [f,df,ddf] = smooth_abs(x,k)
tkx = tanh(k*x);
f   = x.*tkx;
df  = tkx - k*x.*(tkx.^2 - 1);
ddf = 2*k*(k*x.*tkx - 1).*(tkx.^2 - 1);
end

% function [f,df,ddf] = smooth_abs(x,k)
% % ekx = exp(k*x);
% % f   = (2/k)*(log(1+ekx) - log(2)) - x;
% % df  = 2*(ekx./(1+ekx)) - 1;
% % ddf = 2*k*(ekx./(1+ekx).^2);
% 
% % f   = abs(x) + (2/k)*log((1+exp(-k*abs(x)))/2);
% emax = 700; % log(realmax("double"))
% ekx = exp(-min(k*abs(x),emax));
% f   = abs(x) + (2/k)*log((1+ekx)/2);
% mask = x>=0;
% df = zeros(size(x));
% df(mask)  = 2./(1+ekx(mask)) - 1;
% % ekxp = exp(-max(k*x(mask),-emax));
% df(~mask) = 2*ekx(~mask)./(1+ekx(~mask)) - 1;
% ddf = 2*k*df.*(1-df);
% end