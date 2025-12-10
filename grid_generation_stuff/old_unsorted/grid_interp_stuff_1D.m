%% 1D geomspace inverse problem (04/04/2023)
clc; clear; close all;

r0 = 0.5;
rinf = 20;
alpha = 1.1648336;

Nfine = [129,33];

x0 = linspace(r0,rinf,Nfine(2));

r = geomspace(r0,rinf,alpha,Nfine(2));

phi0 = @(x) interp1( x0, r , x, "pchip" );
phi1 = @(N) interp1( x0, r , linspace(r0,rinf,N), "pchip" );

% r2 = phi(linspace(r0,rinf,Nfine(2)));
Nlevels = 9;

N = 2.^(1:Nlevels) + 1;

R = cell(Nlevels,1);
for i = 1:Nlevels
    R{i} = phi1(N(i));
end

hold on
% for i = 1:Nlevels
%     plot(linspace(0,1,N(i)-2), diff(R{i}(2:end))./diff(R{i}(1:end-1)) )
% end

for i = 1:Nlevels
    plot(linspace(0,1,N(i)), R{i} )
end
plot(linspace(0,1,N(5)), R{5},'g:' )