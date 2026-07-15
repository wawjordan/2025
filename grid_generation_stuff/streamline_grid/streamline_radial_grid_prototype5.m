%% streamline/radial_grid prototype 5 (07/14/2026)
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
airfoil_inputs.kappa   = 0.1;
airfoil_inputs.tau     = 0.0;
airfoil_inputs.vinf    = 75.0;
airfoil_inputs.rhoinf  = 1.0;
airfoil_inputs.pinf    = 100000.0;
airfoil_inputs.gamma   = 1.4;
airfoil_inputs.alpha   = 0; % (degrees)
airfoil_inputs.rho_ref = 1.0;
airfoil_inputs.p_ref   = 100000.0;


boundary_distance = 500;
n_body_pts        = 129;
n_wake_pts        = 65;
n_transition      = 1;
jmax              = 129;
n_layers_direct   = 5;
% n_layers_ode      = jmax-n_layers_direct;
n_layers_ode      = 15;
spline_order      = 5;

AR_LE             = 1.0;
AR_TE             = 1.0;
delta_LE          = 0.01;
airfoil = local_airfoil_generator(airfoil_inputs);
[ theta_u,theta_l, psi_1, psi_FF, phi_TE, h_LE, h_TE, delta_TE ] = get_streamline_geometry( ...
                                                            airfoil,  ...
                                                            delta_LE, ...
                                                            AR_LE,    ...
                                                            AR_TE,    ...
                                                            boundary_distance );
[ x_airfoil, y_airfoil, theta_airfoil, phi_airfoil, F1, dF1, ddF1 ] = get_airfoil_coordinates( ...
                                                            airfoil,    ...
                                                            delta_LE,   ...
                                                            delta_TE,   ...
                                                            n_body_pts, ...
                                                            n_transition );


if n_wake_pts>1
    [x_wake,y_wake,~] = get_wake_pts(airfoil,theta_airfoil,n_wake_pts,boundary_distance);
    imax = 2*(n_wake_pts-1) + n_body_pts;
    itmp1 = n_wake_pts;
    itmp2 = n_wake_pts+n_body_pts-1;
    x_all = [flip(x_wake(2:end));x_airfoil;x_wake(2:end)];
    y_all = [flip(y_wake(2:end));y_airfoil;y_wake(2:end)];
else
    imax = n_body_pts;
    itmp1 = 1;
    itmp2 = imax;
    x_all = x_airfoil;
    y_all = y_airfoil;
end

hold on
plot(x_all,y_all,'r.-')
axis equal

% hold on
% plot(x_airfoil,y_airfoil,'r.-')
% axis equal

% allocate
X = zeros(imax,jmax);
Y = zeros(imax,jmax);

X(:,1) = x_all;
Y(:,1) = y_all;

% geometric spacing in psi
psi_vec = my_geomspace(jmax,0,xmax=psi_FF,dx0=psi_1);

% determine indices for transitioning between streamlines and radial
[~,idx_l] = min(abs(theta_airfoil-theta_l));
[~,idx_u] = min(abs(theta_airfoil-theta_u));

i0 = 1;
i1 = 1;
i2 = (itmp1-1) + idx_l;
i3 = (itmp1-1) + idx_u;
i4 = imax;
i5 = imax;

if (n_wake_pts>1)
    i1 = n_wake_pts;
    i4 = imax+1-n_wake_pts;
end

n_middle = idx_u - idx_l + 1;

%% lower surface
for j = 2:1+n_layers_direct
    phi_tmp =  phi_airfoil(1:idx_l);
    psi_tmp = -psi_vec(j);
    z_tmp   = airfoil.z_from_phi_psi(phi_tmp,psi_tmp);
    X(i1:i2,j) = real(z_tmp); Y(i1:i2,j) = imag(z_tmp);
end

xtmp = X(i1:i2,1:1+n_layers_direct);
ytmp = Y(i1:i2,1:1+n_layers_direct);
plot(xtmp,ytmp,'k')
plot(xtmp.',ytmp.','k')
xlim([-1,2])
if ( n_layers_ode>0)
    j1 = 1+n_layers_direct;
    j2 =   n_layers_direct+n_layers_ode;
    psi_tmp = -psi_vec(j1+1:j2);
    [X(i1:i2,j1+1:j2),Y(i1:i2,j1+1:j2)] = ode_phi_extrude_mesh( airfoil, X(i1:i2,j1),Y(i1:i2,j1),psi_tmp,spline_order);
    xtmp = X(i1:i2,j1:j2);
    ytmp = Y(i1:i2,j1:j2);
    plot(xtmp,ytmp,'b')
    plot(xtmp.',ytmp.','b')
end

%% upper surface
for j = 2:1+n_layers_direct
    phi_tmp =  phi_airfoil(idx_u:end);
    psi_tmp = +psi_vec(j);
    z_tmp   = airfoil.z_from_phi_psi(phi_tmp,psi_tmp);
    X(i3:i4,j) = real(z_tmp); Y(i3:i4,j) = imag(z_tmp);
end
xtmp = X(i3:i4,1:1+n_layers_direct);
ytmp = Y(i3:i4,1:1+n_layers_direct);
plot(xtmp,ytmp,'k')
plot(xtmp.',ytmp.','k')

if ( n_layers_ode>0)
    j1 = 1+n_layers_direct;
    j2 =   n_layers_direct+n_layers_ode;
    psi_tmp = +psi_vec(j1+1:j2);
    [X(i3:i4,j1+1:j2),Y(i3:i4,j1+1:j2)] = ode_phi_extrude_mesh( airfoil, X(i3:i4,j1),Y(i3:i4,j1),psi_tmp,spline_order);
    xtmp = X(i3:i4,j1:j2);
    ytmp = Y(i3:i4,j1:j2);
    plot(xtmp,ytmp,'b')
    plot(xtmp.',ytmp.','b')
end

% for j = 2:j2
%     dir = -1;
%     x0 = X(i3,j);
%     y0 = Y(i3,j);
%     phim1 = airfoil.phi_from_xy(X(i3+1,j),y0);
%     phi0 = airfoil.phi_from_xy(x0,y0);
%     phi1 = airfoil.phi_from_xy(x0-1,y0);
%     n_phi = ceil(abs(phi1-phi0)/abs(phi0-phim1));
%     phivec = linspace(phi0,phi1,n_phi);
%     [x,y,xsp,ysp] = ode_psi_extrude_streamline(airfoil,x0,y0,phivec,dir,spline_order);
%     plot(x,y)
%     plot(x0,y0,'rx')
% 
%     x0 = X(i2,j);
%     y0 = Y(i2,j);
%     phim1 = airfoil.phi_from_xy(X(i2-1,j),y0);
%     phi0 = airfoil.phi_from_xy(x0,y0);
%     phi1 = airfoil.phi_from_xy(x0-1,y0);
%     n_phi = ceil(abs(phi1-phi0)/abs(phi0-phim1));
%     phivec = linspace(phi0,phi1,n_phi);
%     [x,y,xsp,ysp] = ode_psi_extrude_streamline(airfoil,x0,y0,phivec,dir,spline_order);
%     plot(x,y)
% end

for j = 2:j2
    dir = -1;
    x0 = X(i3,j);
    y0 = Y(i3,j);
    xm1 = X(i3+1,j);
    x1 = 0;
    psi = airfoil.psi_from_xy(x0,y0);
    n_x = ceil(abs(x1-x0)/abs(x0-xm1));
    x = linspace(x0,x1,n_x);
    y = airfoil.psi_y_from_x(psi,x,y0);
    h = plot(x,y);

    z_l = x(end) + 1i*y(end);
    z_lm1 = x(end-1) + 1i*y(end-1);
    [th_l,r_l] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z_l+airfoil.zLE));
    [th_lm1,~] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z_lm1+airfoil.zLE));
    

    x0 = X(i2,j);
    y0 = Y(i2,j);
    xm1 = X(i2-1,j);
    x1 = 0;
    psi = airfoil.psi_from_xy(x0,y0);
    n_x = ceil(abs(x1-x0)/abs(x0-xm1));
    x = linspace(x0,x1,n_x);
    y = airfoil.psi_y_from_x(psi,x,y0);
    plot(x,y,'Color',h.Color)

    z_u = x(end) + 1i*y(end);
    [th_u,r_u] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z_u+airfoil.zLE));
    th_u = mod(th_u,2*pi);
    th_l = mod(th_l,2*pi);

    th_mid = mod(airfoil.thetaCmax,2*pi);
    z_mid = 0.5*(r_l+r_u)*airfoil.unit_normal_cmplx(th_mid);
    [th_mid2,r_mid] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z_mid+airfoil.zLE));

    n_th = 33;
    th = linspace(th_u,th_l,n_th);
    % r  = linspace(r_u,r_l,n_th);
    thv = [th_l,th_mid,th_u];
    rv  = [r_l,r_mid,r_u];
    r = spline(thv,rv,th);
    zth = airfoil.z_from_theta_r(th,r);
    x = real(zth); y = imag(zth);
    plot(x,y,'Color',h.Color)
end


% for j = 2:n_layers_direct+n_layers_ode
%     h1 = sqrt( ( X(i2,j) - X(i2,j-1) )^2 + ( Y(i2,j) - Y(i2,j-1) )^2 );
%     h2 = sqrt( ( X(i3,j) - X(i3,j-1) )^2 + ( Y(i3,j) - Y(i3,j-1) )^2 );
%     h = linspace(h1,h2,n_middle);
%     % xtmp = X(i1:i2,j-1);
%     % ytmp = Y(i1:i2,j-1);
%     [X(i2:i3,j),Y(i2:i3,j)] = extrude_surface_pts(X(i2:i3,j-1),Y(i2:i3,j-1),1,n_middle,h);
% end
% j1 = 1+n_layers_direct;
% j2 =   n_layers_direct+n_layers_ode;
% xtmp = X(i2:i3,1:j2);
% ytmp = Y(i2:i3,1:j2);
% plot(xtmp,ytmp,'m')
% plot(xtmp.',ytmp.','m')


%% leading edge
% for j = 2:n_layers_direct+n_layers_ode
%     h1 = sqrt( ( X(i2,j) - X(i2,j-1) )^2 + ( Y(i2,j) - Y(i2,j-1) )^2 );
%     h2 = sqrt( ( X(i3,j) - X(i3,j-1) )^2 + ( Y(i3,j) - Y(i3,j-1) )^2 );
%     h = linspace(h1,h2,n_middle);
%     % xtmp = X(i1:i2,j-1);
%     % ytmp = Y(i1:i2,j-1);
%     [X(i2:i3,j),Y(i2:i3,j)] = extrude_surface_pts(X(i2:i3,j-1),Y(i2:i3,j-1),1,n_middle,h);
% end
% j1 = 1+n_layers_direct;
% j2 =   n_layers_direct+n_layers_ode;
% xtmp = X(i2:i3,1:j2);
% ytmp = Y(i2:i3,1:j2);
% plot(xtmp,ytmp,'m')
% plot(xtmp.',ytmp.','m')


if (n_wake_pts>1)
    j1 = 1;
    j2 = n_layers_direct+n_layers_ode;

    psi_tmp = -psi_vec(j1+1:j2);
    i1 = i1-1;
    [X(i0:i1,j1+1:j2),Y(i0:i1,j1+1:j2)] = ode_phi_extrude_mesh( airfoil, X(i0:i1,j1),Y(i0:i1,j1),psi_tmp,spline_order);
    i1 = i1+1;

    j1 = j1+1;
    psi_tmp = -psi_vec(j1+1:j2);
    [X(i1:i1,j1+1:j2),Y(i1:i1,j1+1:j2)] = ode_phi_extrude_mesh( airfoil, X(i1:i1,j1),Y(i1:i1,j1),psi_tmp,spline_order);
    j1 = j1-1;
    xtmp = X(i0:i1,j1:j2);
    ytmp = Y(i0:i1,j1:j2);
    
    plot(xtmp,ytmp,'b')
    plot(xtmp.',ytmp.','b')


    psi_tmp = +psi_vec(j1+1:j2);
    i4 = i4+1;
    [X(i4:i5,j1+1:j2),Y(i4:i5,j1+1:j2)] = ode_phi_extrude_mesh( airfoil, X(i4:i5,j1),Y(i4:i5,j1),psi_tmp,spline_order);
    i4 = i4-1;

    j1 = j1+1;
    psi_tmp = +psi_vec(j1+1:j2);
    [X(i4:i4,j1+1:j2),Y(i4:i4,j1+1:j2)] = ode_phi_extrude_mesh( airfoil, X(i4:i4,j1),Y(i4:i4,j1),psi_tmp,spline_order);
    j1 = j1-1;

    xtmp = X(i4:i5,j1:j2);
    ytmp = Y(i4:i5,j1:j2);
    plot(xtmp,ytmp,'b')
    plot(xtmp.',ytmp.','b')
end


figure;
xtmp = X(:,j1:j2);
ytmp = Y(:,j1:j2);
hold on
plot(xtmp,ytmp,'k')
plot(xtmp.',ytmp.','k')
plot(x_airfoil,y_airfoil,'r.-')
axis equal
xlim([-1,2])
ylim([-1.5,1.5])

j1;

function [theta_u,theta_l,psi_1,psi_FF,phi_TE,h_LE,h_TE,delta_TE] = get_streamline_geometry( airfoil, delta_LE, AR_LE, AR_TE, boundary_distance )
% given streamwise spacing at LE, target AR (streamwise/normal)
% for LE and TE, need to find trailing edge spacing

% approximate cell height at LE
h_LE = delta_LE/AR_LE;

% compute the angles where phi=0

% upper surface
theta_u = fzero( @(theta) real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta) ) ), [0,pi] );
z_u = airfoil.z_from_zeta(airfoil.zeta_from_theta(theta_u));
z_u = (z_u-airfoil.zLE)/airfoil.chord;
z_u = z_u + h_LE*airfoil.unit_normal_cmplx(theta_u);
% find the corresponding value of psi
psi_u = airfoil.psi_from_xy( real(z_u), imag(z_u) );

% lower surface
theta_l = fzero( @(theta) real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta) ) ), [pi,2*pi] );
z_l = airfoil.z_from_zeta(airfoil.zeta_from_theta(theta_l));
z_l = (z_l-airfoil.zLE)/airfoil.chord;
z_l = z_l + h_LE*airfoil.unit_normal_cmplx(theta_l);
% find the corresponding value of psi
psi_l = airfoil.psi_from_xy( real(z_l), imag(z_l) );

% let's average them to maybe get the height at LE closer to target ...
psi_1 = sign(psi_u)*0.5*( abs(psi_u) + abs(psi_l) );

% from here to the leading edge, we are assuming approximately constant cell height

% compute the height at the trailing edge for this psi
zTE = airfoil.zTE;
zTE = (zTE-airfoil.zLE)/airfoil.chord;
phi_TE = airfoil.phi_from_xy(real(zTE),imag(zTE));
z3 = airfoil.z_from_phi_psi(phi_TE,psi_1);
h_TE = abs(z3-zTE);

% compute the corresponding surface spacing
delta_TE = AR_TE*h_TE;

% compute outer boundary streamline

% construct a circle, centered at the trailing edge, and find a location
% such that it has the same potential (phi) as at the trailing edge
%       The termination of this potential line will then be exactly
%       boundary distance away from the trailing edge
%       (at least at that terminal point)

xfun = @(theta) real(zTE) + boundary_distance.*cos(theta);
yfun = @(theta) imag(zTE) + boundary_distance.*sin(theta);
objfun = @(theta) phi_TE-airfoil.phi_from_xy(xfun(theta),yfun(theta));
theta = fzero(@(theta)objfun(theta),pi/2);
x_ff = xfun(theta);
y_ff = yfun(theta);
psi_FF = airfoil.psi_from_xy(x_ff,y_ff);
end

function [x,y,theta,phi,F1,dF1,ddF1] = get_airfoil_coordinates( airfoil, delta_LE, delta_TE, n_pts, n_transition )
% get body coordinates, angle in zeta plane, and stretching function

% get parametric location of maximum curvature on the airfoil
% which will approximately be the leading edge
tLE = airfoil.thetaCmax/(2*pi);

% arc length
L = airfoil.airfoil_arc_length(0,1);

% get normalized arc length to point of max curvature
t0 = airfoil.airfoil_arc_length(0,tLE)/L;

% nondimensionalize spacings
cxL   = airfoil.chord/L;

dLE = delta_LE*cxL;
dTE = delta_TE*cxL;

% get buffer distance over which to blend the two functions
offset = ( (n_transition-1)/(n_pts-1) );

% generate asymmetric stretching function for airfoil surface
[F1,dF1,ddF1] = hermite_blend_2_vinokur_asym(n_pts,t0,dLE,dTE,offset,true);
[x,y,theta] = airfoil.output_airfoil_coords1(n_pts,F1);

% flip to match the clockwise convention
x = flip(x);
y = flip(y);
theta = flip(theta);

% potential along airfoil
phi   = real(airfoil.F_cylinder(airfoil.zeta_from_theta(theta)));
end

function [x,y,phi_ff] = get_wake_pts(airfoil,theta_airfoil,n_wake_pts,boundary_distance)
% Set an exponential stretching rate for the wake region

phi0_l = real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta_airfoil(1)) ) );
phi1_l = real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta_airfoil(2)) ) );

phi0_u = real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta_airfoil(end)) ) );
phi1_u = real( airfoil.F_cylinder( airfoil.zeta_from_theta(theta_airfoil(end-1)) ) );

delta_phi = sign(phi0_u)*0.5*( abs(phi1_l-phi0_l) + abs(phi1_u-phi0_u) );

phi0 = sign(phi0_u)*0.5*( abs(phi0_l)+abs(phi0_u) );

% find phi at outflow boundary
zTE = (airfoil.zTE-airfoil.zLE)/airfoil.chord;
xfun = @(theta) real(zTE) + boundary_distance.*cos(theta);
yfun = @(theta) imag(zTE) + boundary_distance.*sin(theta);
objfun = @(theta) airfoil.psi_from_xy(xfun(theta),yfun(theta));
theta = fzero(@(theta)objfun(theta),0.1);
x_ff = xfun(theta);
y_ff = yfun(theta);
phi_ff = airfoil.phi_from_xy(x_ff,y_ff);

%   do exponential stretching in phi along the stagnation streamline
[phi_vec,~,~] = my_geomspace( n_wake_pts, phi0, xmax=phi_ff, dx0=delta_phi );
phi_vec = phi_vec(:);
options = optimset('TolFun',1e-15,'TolX',1e-17);
x_guess = @(phi) real(zTE) + boundary_distance*(phi-phi0)./(phi_ff-phi0);
% solve for x and y
x = arrayfun( @(phi) fzero( @(x) phi - airfoil.phi_from_xy(x,airfoil.psi_y_from_x(0,x,0)), x_guess(phi), options ),phi_vec);
y = airfoil.psi_y_from_x(0,x,0);
end

function [x,y] = ode_phi_extrude_mesh(airfoil,xvec,yvec,psivec,spline_order)
% numerically integrate equipotential lines and interpolate results with spline
ni = numel(xvec);
nj = numel(psivec);
x = zeros(ni,nj);
y = zeros(ni,nj);
for i = 1:ni
    x0 = xvec(i);
    y0 = yvec(i);
    % factor of 1.1 to make sure we have enough range to interpolate
    [x_,y_,~,ys,~,~,~] = ode_phi_line(airfoil,x0,y0,1.1*psivec(end));
    zeta = ys(:,1) + 1i*ys(:,2);
    psi_ode = imag( airfoil.F_cylinder(zeta) );
    xsp = spapi(spline_order,psi_ode,x_);
    x(i,:) = fnval(xsp,psivec);
    ysp = spapi(spline_order,psi_ode,y_);
    y(i,:) = fnval(ysp,psivec);
end
end

function [x,y,xsp,ysp] = ode_psi_extrude_streamline(airfoil,x,y,phivec,dir,spline_order)
% numerically integrate streamlines and interpolate results with spline
% factor of 1.1 to make sure we have enough range to interpolate
[x_,y_,ts,ys,te,ye,ie] = ode_psi_line(airfoil,x,y,1.1*phivec(end),dir);
zeta    = ys(:,1) + 1i*ys(:,2);
phi_ode = real( airfoil.F_cylinder(zeta) );
xsp = spapi(spline_order,phi_ode,x_);
x   = fnval(xsp,phivec);
ysp = spapi(spline_order,phi_ode,y_);
y   = fnval(ysp,phivec);
end

function [x,y,ts,ys,te,ye,ie] = ode_phi_line(airfoil,x0,y0,psi1)
% numerically integrate equipotential line
% allows for adjusting the integration time automatically
z0 = x0 + 1i*y0;
z0 = z0*airfoil.chord + airfoil.zLE;
zeta0 = airfoil.zeta_from_z(z0);
phi0 = real( airfoil.F_cylinder(zeta0) );
z1 = airfoil.z_from_phi_psi(phi0,psi1);
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

function [x,y,ts,ys,te,ye,ie] = ode_psi_line(airfoil,x0,y0,phi1,dir)
% numerically integrate equipotential line
% allows for adjusting the integration time automatically
z0 = x0 + 1i*y0;
z0 = z0*airfoil.chord + airfoil.zLE;
zeta0 = airfoil.zeta_from_z(z0);
psi0 = imag( airfoil.F_cylinder(zeta0) );
z1 = airfoil.z_from_phi_psi(phi1,psi0);
zeta1 = airfoil.zeta_from_z(z1);
vmag0 = abs(airfoil.w_airfoil(zeta0));
vmag1 = abs(airfoil.w_airfoil(zeta1));

% very rough guess for the integration time
tstop = 2*abs(z1-z0)/abs(vmag1-vmag0);
% s     = sign(psi1);
s = dir;
tspan = s*[0,tstop];
ys0 = [real(zeta0), imag(zeta0)];
[ts,ys,te,ye,ie] = ode_psi_line_step(airfoil,ys0,tspan,phi1);

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
        [ts,ys,te,ye,ie] = ode_psi_line_step(airfoil,ys0,tspan,phi1);
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

    function dydt = phi_integrate(~,y,airfoil)
        % dydt = psi_integrate(0,y,airfoil);
        % tmp = dydt(1);
        % dydt(1) = -dydt(2);
        % dydt(2) = tmp;
        zeta = y(1) + 1i*y(2);
        w = airfoil.w_cylinder(zeta);
        dydt = zeros(2,1);
        dydt(1) = imag(w);
        dydt(2) = real(w);
    end

    function [position,isterminal,direction] = phiEventFcn(~,y,airfoil,psi_target)
      z = airfoil.z_from_zeta(y(1)+1i*y(2));
      z = (z - airfoil.zLE)/airfoil.chord;
      position   = airfoil.psi_from_xy(real(z),imag(z))-psi_target; % The value that we want to be zero
      isterminal = 1;  % Halt integration 
      direction  = 0;   % The zero can be approached from either direction
    end
end

function [ts,ys,te,ye,ie] = ode_psi_line_step(airfoil,ys0,tspan,phi_stop)
    options = odeset('RelTol',1e-10,'AbsTol',1e-12,'NormControl','on');
    options = odeset(options,'Events',@(t,y)psiEventFcn(t,y,airfoil,phi_stop));
    [ts,ys,te,ye,ie] = ode23s(@(t,y) psi_integrate(t,y,airfoil), tspan, ys0, options);

    function dydt = psi_integrate(~,y,airfoil)
        zeta = y(1) + 1i*y(2);
        w = airfoil.w_cylinder(zeta);
        dydt = zeros(2,1);
        dydt(1) = real(w);
        dydt(2) =-imag(w);
    end

    function [position,isterminal,direction] = psiEventFcn(~,y,airfoil,phi_target)
      z = airfoil.z_from_zeta(y(1)+1i*y(2));
      z = (z - airfoil.zLE)/airfoil.chord;
      position   = airfoil.phi_from_xy(real(z),imag(z))-phi_target; % The value that we want to be zero
      isterminal = 1;  % Halt integration 
      direction  = 0;   % The zero can be approached from either direction
    end
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