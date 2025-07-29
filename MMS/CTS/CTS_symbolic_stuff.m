%% CTS Manufactured Solution (revisited) 07/29/2025
clc; clear; close all;
b      = sym('b',[22 1]); % constant coefficients
coeffs = sym('a',[22 5]); % constant coefficients
syms l c mu gamma PrL real
syms x y z t real
f(x,y,z,t) = b( 1) ...
           + b( 2)*sin( pi * ( b( 3)*x   / l   + b( 4) ) ) ...
           + b( 5)*sin( pi * ( b( 6)*y   / l   + b( 7) ) ) ...
           + b( 8)*sin( pi * ( b( 9)*z   / l   + b(10) ) ) ...
           + b(11)*sin( pi * ( b(12)*x*y / l^2 + b(13) ) ) ...
           + b(14)*sin( pi * ( b(15)*y*z / l^2 + b(16) ) ) ...
           + b(17)*sin( pi * ( b(18)*x*z / l^2 + b(19) ) ) ...
           + b(20)*sin( pi * ( b(21)*t*c / l   + b(22) ) );
rho = subs(f,b,coeffs(:,1));
xi  = {x; y; z};
u   = { subs(f,b,coeffs(:,2)); ...
        subs(f,b,coeffs(:,3)); ...
        subs(f,b,coeffs(:,4)) };
p   = subs(f,b,coeffs(:,5));

S = cell(3,3);
for j = 1:3
    for i = 1:3
        S{i,j} = (1/2) * ( diff(u{i},xi(j))  + diff(u{j},xi(i)));
    end
end

tau = cellfun(@(x) 2*mu*x, S, 'UniformOutput',false);

eqn = cell(5,1);

eqn{1} = diff(rho,t);
for i = 1:3
    eqn{1} = eqn{1} + diff(rho*u{i},xi{i});
end

for i = 1:3
    eqn{1+i} = diff(rho*u{i},t) + diff(p,xi{i});
    for j = 1:3
        eqn{1+i} = eqn{1+i} + diff(rho*u{j}*u{i},xi{j}) - diff(tau{j,i},xi{j});
    end
end

rhoE = p/(gamma-1) + (1/2) * rho * (u{1}^2 + u{2}^2 + u{3}^2);
rhoH = p/(gamma-1) + p + (1/2) * rho * (u{1}^2 + u{2}^2 + u{3}^2);
eqn{5} = diff(rhoE,t);

for j = 1:3
    tmp = 0;
    for i = 1:3
        tmp = tmp + u{i}*tau{i,j};
    end
    eqn{5} = eqn{5} + diff(u{j}*rhoH,xi{j}) - diff(tmp,xi{j}) - mu/PrL*diff(rhoH/rho,xi{j});
end

matlabFunction(rho ,"File","dens_mms","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(u{1},"File","uvel_mms","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(u{2},"File","vvel_mms","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(u{3},"File","wvel_mms","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(p   ,"File","pres_mms","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});

matlabFunction(eqn{1},"File","mass_source","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(eqn{2},"File","xmtm_source","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(eqn{3},"File","ymtm_source","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(eqn{4},"File","zmtm_source","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
matlabFunction(eqn{5},"File","ener_source","Vars",{x,y,z,t,coeffs,l,c,mu,gamma,PrL});
