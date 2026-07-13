function [x,y,theta,F1,psi_1,phi_TE,h_LE,h_TE] = get_airfoil_streamline_body_coordinates( airfoil, ...
                                                         n_body_pts, ...
                                                         n_transition, ...
                                                         delta_LE,     ...
                                                         AR_LE, AR_TE )

% given streamwise spacing at LE, target AR (streamwise/normal)
% for LE and TE, need to find trailing edge spacing

% approximate cell height at LE
h_LE = delta_LE/AR_LE;

% get the corresponding point
z1 = h_LE*airfoil.unit_normal_cmplx(airfoil.thetaCmax);

% find radius that this corresponds to in the zeta plane
[~,r1] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z1+airfoil.zLE));

% compute the angle where phi=0
z1 = airfoil.z_from_phi_psi(0,0);
[theta0,~] = airfoil.theta_from_zeta(airfoil.zeta_from_z(airfoil.chord*z1+airfoil.zLE));

% compute the point in the z plane
z2 = airfoil.z_from_theta_r(theta0,r1);

% find the corresponding value of psi
psi_1 = airfoil.psi_from_xy(real(z2),imag(z2));


% r = airfoil.r_from_psi_theta(psi_1,theta0(1));

% compute the height at the trailing edge for this psi
zTE = airfoil.zTE;
zTE = (zTE-airfoil.zLE)/airfoil.chord;

% calculate streamline height at TE
phi_TE = airfoil.phi_from_xy(real(zTE),imag(zTE));
z2 = airfoil.z_from_phi_psi(phi_TE,psi_1);
h_TE = abs(z2-zTE);

% compute the corresponding surface spacing
delta_TE = AR_TE*h_TE;

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
offset = ( (n_transition-1)/(n_body_pts-1) );

% generate asymmetric stretching function for airfoil surface
[F1,~,~] = hermite_blend_2_vinokur_asym(n_body_pts,t0,dLE,dTE,offset,true);
[x,y,theta] = airfoil.output_airfoil_coords1(n_body_pts,F1);

% flip to match the clockwise convention
x = flip(x);
y = flip(y);
theta = flip(theta);
end