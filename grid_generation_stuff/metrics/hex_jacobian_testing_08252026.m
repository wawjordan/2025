%% Hex Jac Testing (08/25/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
grid_file = fullfile(parent_dir,'kt.grd');
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

% x = [-0.304717  1.602827 -0.868254  0.752390  0.574226  0.952229  0.624529  0.590867].';
% y = [ 0.878864  0.866947  0.381305  1.032608  1.391918 -0.079251  1.400126  0.804365].';
% z = [ 0.863041  0.377975 -0.637143  1.120156  0.388743  1.003006  0.135241  1.274976].';
% Jtest = [ 1.865  1.094  0.975  0.712  0.352  0.797  0.284  0.243 ].';
% 
% pts = [x.';y.';z.'];
% 
% J1 = det( hex_jacobian(x,y,z, 0, 0, 0) );
% J2 = det( hex_jacobian(x,y,z, 1, 0, 0) );
% J3 = det( hex_jacobian(x,y,z, 0, 1, 0) );
% J4 = det( hex_jacobian(x,y,z, 1, 1, 0) );
% J5 = det( hex_jacobian(x,y,z, 0, 0, 1) );
% J6 = det( hex_jacobian(x,y,z, 1, 0, 1) );
% J7 = det( hex_jacobian(x,y,z, 0, 1, 1) );
% J8 = det( hex_jacobian(x,y,z, 1, 1, 1) );
% Jvals = [J1; J2; J3; J4; J5; J6; J7; J8];

x = [0,1,0,1,0,1,0,1].'-0.5;
y = [0,0,1,1,0,0,1,1].'-0.5;
z = [0,0,0,0,1,1,1,1].'-0.5;
r1 = vector2vector(rand(3,1),rand(3,1));
% pts = r1*[x.';y.';z.'];
pts = r1*[100*x.';2*y.';z.'];
pts = pts + 0.5*rand(size(pts));
x = pts(1,:).';
y = pts(2,:).';
z = pts(3,:).';

ptmap = [1,2,4,3,1,  5,6,8,7,5,  7,3,4,8,6,2];

pts = [x(ptmap).';y(ptmap).';z(ptmap).'];


M = hex_jacobian(x,y,z, 0.5, 0.5, 0.5);
R1 = get_rotation_matrix(M);
R2 = get_rotation_matrix_2(M);

R1_ = get_rotation_matrix_vv(M,1);
R2_ = get_rotation_matrix_vv(M,2);
R3_ = get_rotation_matrix_vv(M,3);

hold on
plot3(pts(1,:),pts(2,:),pts(3,:),'k')
plot3(pts(1,1),pts(2,1),pts(3,1),'ko')
% pts2 = R1.'*pts;
% plot3(pts2(1,:),pts2(2,:),pts2(3,:),'r')
% plot3(pts2(1,1),pts2(2,1),pts2(3,1),'ro')

pts2 = M\pts;
plot3(pts2(1,:),pts2(2,:),pts2(3,:),'r')
plot3(pts2(1,1),pts2(2,1),pts2(3,1),'ro')

% pts3 = R1_.'*pts;
pts3 = R1.'*pts;
plot3(pts3(1,1),pts3(2,1),pts3(3,1),'go')
plot3(pts3(1,:),pts3(2,:),pts3(3,:),'g')
axis equal


function val = vmag(v)
val = sqrt( sum(v.^2) );
end

function l = get_jac_length(M,dir)
l = vmag(M(:,dir));
end

function zeta = get_jac_size(M)
l1 = get_jac_length(M,1);
l2 = get_jac_length(M,2);
l3 = get_jac_length(M,3);
zeta = l1*l2*l3;
end

function v = get_jac_vol(M)
v = det(M);
end

function vec = get_jac_unit_vec(M,dir)
vec = M(:,dir)/get_jac_length(M,dir);
end

function phi = get_skew_angle_c(M,dir1,dir2)
cosphi = dot( get_jac_unit_vec(M,dir1), get_jac_unit_vec(M,dir2) );
phi = acos( cosphi );
end

function phi = get_skew_angle_s(M,dir1,dir2)
sinphi = vmag( cross( get_jac_unit_vec(M,dir1), get_jac_unit_vec(M,dir2) ) );
phi = asin( sinphi );
end

function chi = get_dihedral_angle_c0(M)
cr1 = cross( get_jac_unit_vec(M,1), get_jac_unit_vec(M,2) );
cr2 = cross( get_jac_unit_vec(M,1), get_jac_unit_vec(M,3) );
coschi = dot( cr1/vmag(cr1), cr2/vmag(cr2) );
chi = acos(coschi);
end

function chi = get_dihedral_angle_c(M)
phi12 = get_skew_angle_c(M,1,2);
phi13 = get_skew_angle_c(M,1,3);
phi23 = get_skew_angle_c(M,2,3);
coschi = ( cos(phi23) - cos(phi12)*cos(phi13) )/(sin(phi12)*sin(phi13));
chi = acos(coschi);
end

function chi = get_dihedral_angle_s(M)
phi12 = get_skew_angle_c(M,1,2);
phi13 = get_skew_angle_c(M,1,3);
zeta  = get_jac_size(M);
sinchi = det(M)/( zeta*sin(phi12)*sin(phi13) );
chi = asin(sinchi);
end

function rho = get_aspect_ratio(M,dir)
ind = 1:3;
l = zeros(3,1);
l(1) = get_jac_length(M,1);
l(2) = get_jac_length(M,2);
l(3) = get_jac_length(M,3);
rho = l(dir)/sqrt( prod( l(ind~=dir) ) );
end

function theta = get_theta_orientation_c(M)
costheta = M(1,1)/sqrt( M(1,1)^2 + M(2,1)^2 );
theta    = acos(costheta);
end

function theta = get_theta_orientation_s(M)
sintheta = M(2,1)/sqrt( M(1,1)^2 + M(2,1)^2 );
theta    = asin(sintheta);
end

function psi = get_psi_orientation_c(M)
l1 = get_jac_length(M,1);
cospsi = M(3,1)/l1;
psi = acos(cospsi);
end

function psi = get_psi_orientation_s(M)
l1 = get_jac_length(M,1);
sinpsi = sqrt( M(1,1)^2 + M(2,1)^2 )/l1;
psi = acos(sinpsi);
end

function beta = get_beta(M)
x1 = get_jac_unit_vec(M,1);
x2 = get_jac_unit_vec(M,2);
crx12 = cross( x1, x2 );
cosphi12 = dot( x1, x2 );
sinphi12 = vmag( crx12 );
r2 = (x2 - cosphi12*x1)/sinphi12;
r3 = crx12/sinphi12;

sintheta = M(2,1)/sqrt( M(1,1)^2 + M(2,1)^2 );
costheta = M(1,1)/sqrt( M(1,1)^2 + M(2,1)^2 );

aperp = [-sintheta; +costheta; 0];

cosbeta = dot(r2,aperp);
sinbeta = dot(r3,aperp);

beta = acos( cosbeta );
beta = asin( sinbeta );
end

function R = get_rotation_matrix(M)
R = zeros(3,3);
x1 = get_jac_unit_vec(M,1);
x2 = get_jac_unit_vec(M,2);
crx12 = cross( x1, x2 );
cosphi12 = dot( x1, x2 );
sinphi12 = vmag( crx12 );
R(:,1) = x1;
R(:,2) = ( x2 - cosphi12*x1 )/sinphi12;
R(:,3) = crx12/sinphi12;
end

function R = get_rotation_matrix_2(M)
R = zeros(3,3);
ct = M(1,1)/sqrt( M(1,1)^2 + M(2,1)^2 );
st = M(2,1)/sqrt( M(1,1)^2 + M(2,1)^2 );

l1  = vmag(M(:,1));
cp = M(3,1)/l1;
sp = sqrt( M(1,1)^2 + M(2,1)^2 )/l1;

x1 = get_jac_unit_vec(M,1);
x2 = get_jac_unit_vec(M,2);
crx12 = cross( x1, x2 );
cosphi12 = dot( x1, x2 );
sinphi12 = vmag( crx12 );
r2 = (x2 - cosphi12*x1)/sinphi12;
r3 = crx12/sinphi12;
aperp = [-st; +ct; 0];
cb = dot(r2,aperp);
sb = dot(r3,aperp);

R(1,1) = ct*sp;
R(2,1) = st*sp;
R(3,1) = cp;
R(1,2) = -st*cb + ct*cp*sb;
R(2,2) = ct*cb + st*cp*sb;
R(3,2) = -sp*sb;
R(1,3) = -st*sb -ct*cp*cb;
R(2,3) = ct*sb - st*cp*cb;
R(3,3) = sp*cb;
end

function R = get_rotation_matrix_vv(M,dir)
x = zeros(3,1);
x(dir) = 1;
y = M(:,dir);
R = vector2vector(x,y);
end

function Rxy = vector2vector(x,y)
% https://arxiv.org/pdf/1801.01724
% note: u,v in S^n (unit sphere, so they need to be unit length, v /= -u
x = x(:)/vmag(x);
y = y(:)/vmag(y);
yxt = y*x.';
xyt = x*y.';
xxt = x*x.';
yyt = y*y.';
xyd = dot(x,y);
Rxy = eye(numel(x)) + yxt - xyt + (1/(1+xyd)) * ( xyd*(yxt + xyt) - (xxt + yyt) );
% Kxy = yxt - xyt;
% Rxy = eye(numel(x)) + Kxy + (1/(1+xyd)) * (yxt - xyt)^2;
end