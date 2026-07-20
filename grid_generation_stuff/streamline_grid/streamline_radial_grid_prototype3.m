%% streamline/radial_grid prototype 3 (07/12/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

airfoil_inputs = struct();
airfoil_inputs.epsilon = 0.1;
airfoil_inputs.kappa   = 0.0;
airfoil_inputs.tau     = 0.1;
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = 10; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;
boundary_distance = 500;
n_body_pts = 129;
n_wake_pts = 65;
n_transition = 1;
jmax_in    = 1;
AR_LE      =  1.0;
AR_TE      =  1.0;
delta_LE   = 0.001;
airfoil = local_airfoil_generator(airfoil_inputs);
[x,y,theta0,F1,psi_1,phi_TE,h_LE,h_TE] = get_airfoil_streamline_body_coordinates2( airfoil, ...
                                                         n_body_pts, ...
                                                         n_transition, ...
                                                         delta_LE,     ...
                                                         AR_LE, AR_TE );

dj = 20;
di = n_body_pts;

hold on
plot(x,y,'r.-')
axis equal

% find phi at trailing edge, then find coordinates for that phi at the outer boundary
[phi1,x1,y1,~] = min_phi_on_circle(airfoil,x(1),y(1),boundary_distance,pi/2);

% get the streamline that goes through that point
psi_jmax = airfoil.psi_from_xy(x1,y1);

psi0 = my_geomspace(129,0,xmax=psi_jmax,dx0=psi_1);
% phi0 = airfoil.phi_from_xy(x,y);
phi0 = real(airfoil.F_cylinder(airfoil.zeta_from_theta(theta0)));

% [xs,ys] = ode_streamline_points(airfoil,x0,y0,true,psis);
% xs = [x0;xs];
% ys = [y0;ys];
% plot(xs,ys,'go-')




X = zeros(di,dj);
Y = zeros(di,dj);

X(:,1) = x;
Y(:,1) = y;


[~,idx1] = min(abs(theta0-(airfoil.thetaSP+pi/2)));
[~,idx2] = min(abs(theta0-(airfoil.thetaSP-pi/2)));


%% lower surface
phi =  phi0(1:idx1);
psi = -psi0(2:dj);
[PHI,PSI] = ndgrid(phi,psi);
tic;
Z1 = airfoil.z_from_phi_psi(PHI,PSI);
time1=toc;
X1 = real(Z1); Y1 = imag(Z1);
X(1:idx1,2:dj) = X1; Y(1:idx1,2:dj) = Y1;
plot(X1,Y1,'k')
plot(X1.',Y1.','k')

% [Xlo,Ylo] = ode_phi_extrude_mesh(airfoil,X(1:idx1,dj),Y(1:idx1,dj),-psi0(dj+1:30),5);
% plot(Xlo,Ylo,'b')
% plot(Xlo.',Ylo.','b')


% tic;
% for i = 1:idx1
%     x0 = X(i,dj);
%     y0 = Y(i,dj);
%     [x_,y_,ts,ys,te,ye,ie] = ode_phi_line(airfoil,x0,y0,-psi0(30));
%     plot(x_,y_,'b.-')
% end
% time2=toc;
% disp(time1);
% disp(time2);
%% leading edge
z1 = airfoil.z_from_phi_psi(phi0(idx1),psi);
plot(real(z1),imag(z1),'g')
[theta1,r1] = airfoil.theta_from_zeta(airfoil.zeta_from_zs(z1));
phi =  phi0(idx2:end);
psi =  psi0(2:dj);
z2 = airfoil.z_from_phi_psi(phi0(idx2),psi);
plot(real(z2),imag(z2),'g')
[theta2,r2] = airfoil.theta_from_zeta(airfoil.zeta_from_zs(z2));


theta = theta0(idx1:idx2);
[THETA,R] =ndgrid(theta,r1);
% for j = 1:numel(r1)
%     R(:,j) = linspace(r1(j),r2(j),numel(theta));
%     for i = 1:numel(theta)
%         R(i,j) = airfoil.r_from_psi_theta(psi(j),theta(i));
%     end
% end
Z2 = airfoil.zs_from_theta_r(THETA,R);
X2 = real(Z2); Y2 = imag(Z2);
X(idx1:idx2,2:dj) = X2; Y(idx1:idx2,2:dj) = Y2;
plot(X2,Y2,'b')
plot(X2.',Y2.','b')

%% upper surface
[PHI,PSI] = ndgrid(phi,psi);
Z3 = airfoil.z_from_phi_psi(PHI,PSI);
X3 = real(Z3); Y3 = imag(Z3);
X(idx2+1:end,2:dj) = X3(2:end,:); Y(idx2+1:end,2:dj) = Y3(2:end,:);
% plot(X3,Y3,'k')
% plot(X3.',Y3.','k')

plot(X,Y,'k')
plot(X.',Y.','k')


Xp = zeros(di,dj);
Yp = zeros(di,dj);
Xp(:,1) = x;
Yp(:,1) = y;

for j = 2:dj
    h1 = sqrt( ( X(idx1,j) - X(idx1,j-1) )^2 + ( Y(idx1,j) - Y(idx1,j-1) )^2 );
    h2 = sqrt( ( X(idx2,j) - X(idx2,j-1) )^2 + ( Y(idx2,j) - Y(idx2,j-1) )^2 );
    h = zeros(di,1);
    h(1:idx1-1)  = h1;
    h(idx1:idx2) = linspace(h1,h2,idx2-idx1+1);
    h(idx2+1:end)  = h2;
    [Xp(:,j),Yp(:,j)] = extrude_surface_pts(Xp(:,j-1),Yp(:,j-1),1,di,h);
end


plot(Xp,Yp,'r')
plot(Xp.',Yp.','r')

idx1;

function [phi,x,y,theta] = min_phi_on_circle(airfoil,x0,y0,r,theta_guess)
xfun = @(theta) x0 + r.*cos(theta);
yfun = @(theta) y0 + r.*sin(theta);
phi0 = airfoil.phi_from_xy(x0,y0);
objfun = @(theta) phi0-airfoil.phi_from_xy(xfun(theta),yfun(theta));
theta = fzero(@(theta)objfun(theta),theta_guess);
x = xfun(theta);
y = yfun(theta);
phi = airfoil.phi_from_xy(x,y);
end

function [x,y] = ode_phi_extrude_mesh(airfoil,xvec,yvec,psivec,spline_order)
ni = numel(xvec);
nj = numel(psivec);
x = zeros(ni,nj);
y = zeros(ni,nj);
for i = 1:ni
    x0 = xvec(i);
    y0 = yvec(i);
    [x_,y_,~,ys,~,~,~] = ode_phi_line(airfoil,x0,y0,1.1*psivec(end));
    zeta = ys(:,1) + 1i*ys(:,2);
    psi_ode = imag( airfoil.F_cylinder(zeta) );
    xsp = spapi(spline_order,psi_ode,x_);
    x(i,:) = fnval(xsp,psivec);
    ysp = spapi(spline_order,psi_ode,y_);
    y(i,:) = fnval(ysp,psivec);
end
end

function [x,y,ts,ys,te,ye,ie] = ode_phi_line(airfoil,x0,y0,psi1)
z0 = x0 + 1i*y0;
z0 = z0*airfoil.chord + airfoil.zLE;
zeta0 = airfoil.zeta_from_z(z0);

phi0 = real( airfoil.F_cylinder(zeta0) );

z1 = airfoil.z_from_phi_psi(phi0,psi1);
% z1 = x1 + 1i*y1;
% z1 = z1*airfoil.chord + airfoil.zLE;
zeta1 = airfoil.zeta_from_z(z1);
vmag0 = abs(airfoil.w_airfoil(zeta0));
vmag1 = abs(airfoil.w_airfoil(zeta1));

% very rough guess for the integration time
tstop = 2*abs(z1-z0)/abs(vmag1-vmag0);
s     = sign(psi1);
tspan = s*[0,tstop];
ys0 = [real(zeta0), imag(zeta0)];
[ts,ys,te,ye,ie] = ode_phi_line_step(airfoil,ys0,tspan,psi1);


found = false;
maxiter = 10;
iter = 1;
while (~found)
    if iter>maxiter
        error('failed to converge')
    end  
    if ~(isempty(ie))
        found = true;
    else
        zeta_e = ys(end,1) + 1i*ys(end,2);
        ze = airfoil.z_from_zeta(zeta_e);
        ze = (ze-airfoil.zLE)/airfoil.chord;
        dist_ratio = abs(ze-z0)/abs(z1-z0);
        dist_ratio = min(max(dist_ratio,1.1),3);
        tspan = tspan*dist_ratio;
        [ts,ys,te,ye,ie] = ode_phi_line_step(airfoil,ys0,tspan,psi1);
    end
    iter = iter+1;
end
zeta = ys(:,1) + 1i*ys(:,2);
z = airfoil.z_from_zeta(zeta);
z = (z-airfoil.zLE)/airfoil.chord;
x = real(z);
y = imag(z);
end

function [ts,ys,te,ye,ie] = ode_phi_line_step(airfoil,ys0,tspan,psi_stop)
options = odeset('RelTol',1e-10,'AbsTol',1e-12,'NormControl','on');
options = odeset(options,'Events',@(t,y)phiEventFcn(t,y,airfoil,psi_stop));
[ts,ys,te,ye,ie] = ode23s(@(t,y) phi_integrate(t,y,airfoil), tspan, ys0, options);
end

function [x,y,ts,ys,te,ye,ie] = ode_psi_line_step(airfoil,x0,y0,tspan,phi_stop)
options = odeset('RelTol',1e-12,'AbsTol',1e-16,'InitialStep',1e-12,'NormControl','on');
z0 = x0 + 1i*y0;
z0 = z0*airfoil.chord + airfoil.zLE;
zeta0 = airfoil.zeta_from_z_alt(z0);
ys0 = [real(zeta0), imag(zeta0)];
options = odeset(options,'Events',@(t,y)psiEventFcn(t,y,airfoil,phi_stop));
[ts,ys,te,ye,ie] = ode113(@(t,y) psi_integrate(t,y,airfoil), tspan, ys0, options);
zeta = ys(:,1) + 1i*ys(:,2);
z = airfoil.z_from_zeta(zeta);
z = (z-airfoil.zLE)/airfoil.chord;
x = real(z);
y = imag(z);
end


function [x,y,ts,xe,ye,te] = ode_streamline(airfoil,x0,y0,tspan,potential,end_val)
tol = 1e-12;
dt  = 0.1;
step = 1e-6;
maxiter = 100;
options = odeset('RelTol',1e-12,'AbsTol',1e-16,'InitialStep',1e-12,'NormControl','on');
s = 2*(end_val>=0)-1;
z0 = x0 + 1i*y0;
z0 = z0*airfoil.chord + airfoil.zLE;
zeta0 = airfoil.zeta_from_z_alt(z0);


% attempt to integrate
% [success,ts_tmp,ys_tmp,zeta0] = attempt_integration( airfoil, dt, s*step, tol, maxiter, zeta0, options, potential );
% if ~success
%     [success,ts_tmp,ys_tmp,zeta0] = attempt_integration( airfoil, dt, -s*step, tol, maxiter, zeta0, options, potential );
% end
% if (success)
    ys0 = [real(zeta0), imag(zeta0)];
    if potential
        options = odeset(options,'Events',@(t,y)phiEventFcn(t,y,airfoil,end_val));
        [ts,ys,te,ye,~] = ode113(@(t,y) phi_integrate(t,y,airfoil), s*tspan, ys0, options);
    else
        options = odeset(options,'Events',@(t,y)psiEventFcn(t,y,airfoil,end_val));
        [ts,ys,te,ye,~] = ode113(@(t,y) psi_integrate(t,y,airfoil), s*tspan, ys0, options);
    end
%     ts = [ts_tmp(1:end,1);ts_tmp(end,1)+ts];
%     ys = [ys_tmp(1:end,:);ys];
% else
%     ts = ts_tmp;
%     ys = ys_tmp;
% end
zeta = ys(:,1) + 1i*ys(:,2);
z = airfoil.z_from_zeta(zeta);
z = (z-airfoil.zLE)/airfoil.chord;
x = real(z);
y = imag(z);

if nargout>3
    ze = airfoil.z_from_zeta(ye(:,1) + 1i*ye(:,2));
    ze = (ze-airfoil.zLE)/airfoil.chord;
    xe = real(ze);
    ye = imag(ze);
end
end

function [success,ts,ys,zeta0] = attempt_integration( airfoil, dt, step, tol, maxiter, zeta0, options, potential )
    if (potential)
        fstep = step;
    else
        fstep = 1i*step;
    end
    tspan = [0,dt];
    ys0 = [real(zeta0), imag(zeta0)];
    if potential
        [~,ys_tmp] = ode113(@(t,y) phi_integrate(t,y,airfoil),tspan, ys0, options);
    else
        [~,ys_tmp] = ode113(@(t,y) psi_integrate(t,y,airfoil),tspan, ys0, options);
    end
    success = norm(ys_tmp(end,:)-ys_tmp(1,:))>tol;
    ts = linspace(0,dt*maxiter,maxiter).';
    ys = zeros(maxiter,numel(ys0));
    iter = 1;
    while ~success
        if (iter>maxiter)
            return
        end
        F = airfoil.F_cylinder(zeta0);
        % take a step in psi
        F = F + fstep;
        % find new zeta
        theta = airfoil.theta_from_zeta(zeta0);
        % try marching away from the airfoil
        dzeta = step*exp(1i*(theta-airfoil.beta));
        zeta = airfoil.zeta_from_potential(F,zeta0 + dzeta);
        zeta0=zeta;
        ys0 = [real(zeta0), imag(zeta0)];
        ys(iter,:) = ys0;
        if potential
            [~,ys_tmp] = ode113(@(t,y) phi_integrate(t,y,airfoil),tspan, ys0, options);
        else
            [~,ys_tmp] = ode113(@(t,y) psi_integrate(t,y,airfoil),tspan, ys0, options);
        end
        success = norm(ys_tmp(end,:)-ys_tmp(1,:))>tol;
        iter = iter + 1;
    end
    ts = ts(1:iter,1);
    ys = ys(1:iter,:);
end

function [position,isterminal,direction] = psiEventFcn(~,y,airfoil,phi_target)
  z = airfoil.z_from_zeta(y(1)+1i*y(2));
  z = (z - airfoil.zLE)/airfoil.chord;
  position   = airfoil.phi_from_xy(real(z),imag(z))-phi_target; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction  = 0;   % The zero can be approached from either direction
end

function [position,isterminal,direction] = phiEventFcn(~,y,airfoil,psi_target)
  z = airfoil.z_from_zeta(y(1)+1i*y(2));
  z = (z - airfoil.zLE)/airfoil.chord;
  position   = airfoil.psi_from_xy(real(z),imag(z))-psi_target; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction  = 0;   % The zero can be approached from either direction
end


function dydt = phi_integrate(~,y,airfoil)
    dydt = psi_integrate(0,y,airfoil);
    tmp = dydt(1);
    dydt(1) = -dydt(2);
    dydt(2) = tmp;
end

function dydt = psi_integrate(~,y,airfoil)
zeta = y(1) + 1i*y(2);
% w = airfoil.w_airfoil(zeta);
w = airfoil.w_cylinder(zeta);
dydt = zeros(2,1);
dydt(1) = real(w);
dydt(2) =-imag(w);
end

function airfoil = local_airfoil_generator(in)
% epsilon = 0.1;
% kappa   = 0.0;
% tau     = 0.0;
% vinf    = 75.0;
% rhoinf  = 1.0;
% pinf    = 100000.0;
% gamma   = 1.4;
% alpha   = 0; % (degrees)
% rho_ref = 1.0;
% p_ref   = 100000.0;
a_ref   = sqrt(in.gamma*in.p_ref/in.rho_ref);

% nondimensionalize inputs
vinf   = in.vinf/a_ref;
rhoinf = in.rhoinf/in.rho_ref;
pinf   = in.pinf/(in.rho_ref*a_ref^2);

airfoil        = kt_airfoil( in.epsilon, in.kappa, in.tau );
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(in.alpha);
end