%% Navah-Nadarajah MS-1 source terms 08/03/2026
clc; clear; close all;
syms mu gamma PrL real
syms x y z t real
y_ = y - 0.05*sin(2*pi*x);
dydx = 2*pi*0.05*cos(2*pi*x);
rho0 = 1.0;
p0   = 1.0;
uw   = 1.0;
rho = rho0 + y_.^2;
u   = { uw + y_; dydx.*(uw+y); 0*x+0*y};
p   = p0 + y_.^2;
xi  = {x; y; z};

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

matlabFunction(rho ,"File","dens_mms","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(u{1},"File","uvel_mms","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(u{2},"File","vvel_mms","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(u{3},"File","wvel_mms","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(p   ,"File","pres_mms","Vars",{x,y,z,t,mu,gamma,PrL});

matlabFunction(eqn{1},"File","mass_source","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(eqn{2},"File","xmtm_source","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(eqn{3},"File","ymtm_source","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(eqn{4},"File","zmtm_source","Vars",{x,y,z,t,mu,gamma,PrL});
matlabFunction(eqn{5},"File","ener_source","Vars",{x,y,z,t,mu,gamma,PrL});
