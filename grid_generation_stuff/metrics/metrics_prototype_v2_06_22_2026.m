%% Grid Metric Prototyping v2
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

xAR2 = -1;
while xAR2<=0
    [x,y] = sort_points_counterclockwise(rand([4,1]),rand([4,1]));
    xAR2 = aspect_ratio_2(x,y);
end
% x = [0,1,1,0];
% y = [0,0,1,1];

j1 = jac_ind_1(x,y);

% for j = -1:1
%     for i = -1:1
%         A = quad_jacobian(x,y,i,j);
%         A2 = quad_jacobian2(x,y,(i+1)/2,(j+1)/2);
%         max(abs(A-A2),[],'all')
%     end
% end
A = quad_jacobian2(x,y,0.5,0.5);
[zeta,theta,phi,rho] = parameters_from_jacobian(A);
W = parameterized_2D_matrix(zeta,theta,phi,rho);
[rt_zeta,R,Q,D,S,U] = jacobian_decomp(A);
% max(abs(A- rt_zeta*R*Q*D),[],'all')
c = polygon_centroid(x,y);
% max(abs(quad_centroid(x,y) - polygon_centroid(x,y)),[],'all')
% max(abs(quad_centroid_2(x,y,1) - polygon_centroid(x,y)),[],'all')
% max(abs(quad_centroid_2(x,y,2) - polygon_centroid(x,y)),[],'all')

s1  = skewness_1(x,y);
s2  = skewness_2(x,y);
T   = taper(x,y);
xAR1 = aspect_ratio_1(x,y);
xAR2 = aspect_ratio_2(x,y);
xAR3 = aspect_ratio_3(x,y);
% j2 = jac_ind_1(x,y);

plot([x(:);x(1)],[y(:);y(1)])
hold on
x2 = zeros(4,1);
y2 = zeros(4,1);

% tmp = A*inv(D);
% for i = 1:4
%     x2(i) = tmp
% plot([x(:);x(1)],[y(:);y(1)])
% plot([c(1),c(1)+A(1,1)/vmag(A(:,1))],[c(2),c(2)+A(2,1)/vmag(A(:,1))],'r')
% plot([c(1),c(1)+A(1,2)/vmag(A(:,2))],[c(2),c(2)+A(2,2)/vmag(A(:,2))],'b')
plot([c(1),c(1)+A(1,1)],[c(2),c(2)+A(2,1)],'r')
plot([c(1),c(1)+A(1,2)],[c(2),c(2)+A(2,2)],'b')
axis equal

tmp = inv(D)*inv(Q)*inv(R)*[x.';y.']/rt_zeta;
x2 = tmp(1,:).';
y2 = tmp(2,:).';
A2 = quad_jacobian2(x2,y2,0.5,0.5);
c2 = polygon_centroid(x2,y2);
plot([x2;x2(1)],[y2;y2(1)],'--')
plot([c2(1),c2(1)+A2(1,1)],[c2(2),c2(2)+A2(2,1)],'r--')
plot([c2(1),c2(1)+A2(1,2)],[c2(2),c2(2)+A2(2,2)],'b--')


% tmp = inv(A)*[x.';y.'];
tmp = A\[x.';y.'];
x2 = tmp(1,:).';
y2 = tmp(2,:).';
A2 = quad_jacobian2(x2,y2,0.5,0.5);
c2 = polygon_centroid(x2,y2);
plot([x2;x2(1)],[y2;y2(1)],'--')
plot([c2(1),c2(1)+A2(1,1)],[c2(2),c2(2)+A2(2,1)],'r:')
plot([c2(1),c2(1)+A2(1,2)],[c2(2),c2(2)+A2(2,2)],'b:')

function val = vmag(v)
val = sqrt( sum(v.^2) );
end

function val = skewness_1(x,y)
x1   = ( x(2) - x(1) ) + (x(3) - x(4) );
y1   = ( y(2) - y(1) ) + (y(3) - y(4) );
x2   = ( x(3) - x(2) ) + (x(4) - x(1) );
y2   = ( y(3) - y(2) ) + (y(4) - y(1) );
m1m2 = sqrt( x1^2 + y1^2 ) * sqrt( x2^2 + y2^2 );
% val  = abs( (x1*x2+y1*y2)/m1m2 );
val  = 1-abs( (x1*x2+y1*y2)/m1m2 );
end

function val = skewness_2(x,y)
val = 0;
for k = 1:4
    km1  = mod(k-2,4)+1;
    kp1  = mod(k  ,4)+1;
    e1x = x(km1)-x(k);
    e1y = y(km1)-y(k);
    e2x = x(kp1)-x(k);
    e2y = y(kp1)-y(k);
    m1m2 = sqrt( e1x^2 + e1y^2 ) * sqrt( e2x^2 + e2y^2 );
    val = val + 1/sqrt( 1 - ( (e1x*e2x + e1y*e2y)/m1m2 )^2 );
    % stheta = sqrt( 1 - ( (e1x*e2x + e1y*e2y)/m1m2 )^2 );
    % stheta2 = sin(acos( (e1x*e2x + e1y*e2y)/m1m2 ));
    % stheta-stheta2
end
val = 4/val;
end

function val = taper(x,y)
m1  = sqrt( ( ( x(2) - x(1) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(2) - y(1) ) + (y(3) - y(4) ) )^2 );
m2  = sqrt( ( ( x(3) - x(2) ) + (x(4) - x(1) ) )^2 ...
          + ( ( y(3) - y(2) ) + (y(4) - y(1) ) )^2 );
m12 = sqrt( ( ( x(1) - x(2) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(1) - y(2) ) + (y(3) - y(4) ) )^2 );
% val = m12/min( m1, m2 );
val = 1-m12/min( m1, m2 );
end

function val = aspect_ratio_1(x,y)
m1  = sqrt( ( ( x(2) - x(1) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(2) - y(1) ) + (y(3) - y(4) ) )^2 );
m2  = sqrt( ( ( x(3) - x(2) ) + (x(4) - x(1) ) )^2 ...
          + ( ( y(3) - y(2) ) + (y(4) - y(1) ) )^2 );
val = m1/m2;
% val = max( val, 1/val );
val = min( val, 1/val );
end

function val = aspect_ratio_2(x,y)
hmax = max( sqrt( (x(3)-x(1))^2 + (y(3)-y(1))^2 ), ...
            sqrt( (x(4)-x(2))^2 + (y(4)-y(2))^2 ) );
hmag = 0;
minAx2 = realmax;
for k = 1:4
    kp1  = mod(k  ,4)+1;
    kp3  = mod(k+2,4)+1;
    minAx2 = min( minAx2, ...
                ( x(k)   - x(kp3) )*( y(kp1) - y(kp3)) ...
              - ( y(k)   - y(kp3) )*( x(kp1) - x(kp3) ) );
    h2   = (x(kp1)-x(k))^2 + (y(kp1)-y(k))^2;
    hmag = hmag + h2;
    hmax = max( hmax, sqrt(h2) );
end
% val = 0.25*sqrt(2*hmag)*hmax/minAx2;
val = 1/( 0.25*sqrt(2*hmag)*hmax/minAx2 );
end

function val = aspect_ratio_3(x,y)
x1   = ( x(4) - x(1) );
y1   = ( y(4) - y(1) );
x2   = ( x(2) - x(1) );
y2   = ( y(2) - y(1) );
m1m2 = sqrt( x1^2 + y1^2 ) * sqrt( x2^2 + y2^2 );
alpha = acos( (x1*x2+y1*y2)/m1m2 );
x1   = ( x(4) - x(3) );
y1   = ( y(4) - y(3) );
x2   = ( x(2) - x(3) );
y2   = ( y(2) - y(3) );
m1m2 = sqrt( x1^2 + y1^2 ) * sqrt( x2^2 + y2^2 );
alpha = 0.5*( alpha + acos( (x1*x2+y1*y2)/m1m2 ) );
% val = 1/(1-abs(1-(2*alpha)/pi) );
val = 1-abs(1-(2*alpha)/pi);
end

function val = jac_ind_1(x,y)
den = 0;
for n = 1:4
    A = quad_jacobian_knupp(x,y,n);
    A = A'*A;
    den = den + sum(eig(A))/det(A);
end
val = 8/den;
end

function val = jac_ind_2(x,y)
den = 0;
for j = -1:2:1
    for i = -1:2:1
        A = quad_jacobian(x,y,i,j);
        A = A'*A;
        den = den + sum(eig(A))/det(A);
    end
end
val = 8/den;
end

function J = quad_jacobian_knupp(x,y,k)
J = zeros(2,2);
kp1  = mod(k  ,4)+1;
kp3  = mod(k+2,4)+1;
J(1,1) = x(kp1) - x(k);
J(1,2) = x(kp3) - x(k);
J(2,1) = y(kp1) - y(k);
J(2,2) = y(kp3) - y(k);
end

% function J = quad_jacobian(x,y,xi,eta)
% J = zeros(2,2);
% J(1,1) = x(1) * -(1-eta) ...
%        + x(2) *  (1-eta) ...
%        + x(3) *  (  eta) ...
%        + x(4) * -(  eta);
% J(1,2) = x(1) * -(1-xi) ...
%        + x(2) * -(  xi) ...
%        + x(3) *  (  xi) ...
%        + x(4) *  (1-xi);
% 
% J(2,1) = y(1) * -(1-eta) ...
%        + y(2) *  (1-eta) ...
%        + y(3) *  (  eta) ...
%        + y(4) * -(  eta);
% J(2,2) = y(1) * -(1-xi) ...
%        + y(2) * -(  xi) ...
%        + y(3) *  (  xi) ...
%        + y(4) *  (1-xi);
% J = J*[0,-1;1,0];
% % J(1,1) = ( (x(2)-x(1))*(1-eta) ...
% %          + (x(3)-x(4))*  eta   );
% % 
% % J(2,1) = ( (y(2)-y(1))*(1-eta) ...
% %          + (y(3)-y(4))*  eta   );
% % 
% % J(1,2) = ( (x(4)-x(1))*(1-xi) ...
% %          + (x(3)-x(2))*   xi   );
% % 
% % J(2,2) = ( (y(4)-y(1))*(1-xi) ...
% %          + (y(3)-y(2))*   xi   );
% end

function J = quad_jacobian(x,y,xi,eta)
J = zeros(2,2);
v1 = [x(1),y(1)];
v2 = [x(2),y(2)];
v3 = [x(3),y(3)];
v4 = [x(4),y(4)];
% P1 = 0.25*(1-xi).*(1-eta);
% P2 = 0.25*(1+xi).*(1-eta);
% P3 = 0.25*(1+xi).*(1+eta);
% P4 = 0.25*(1-xi).*(1+eta);
dP1dxi = -0.25*(1-eta);
dP2dxi =  0.25*(1-eta);
dP3dxi =  0.25*(1+eta);
dP4dxi = -0.25*(1+eta);

dP1deta = -0.25*(1-xi);
dP2deta = -0.25*(1+xi);
dP3deta =  0.25*(1+xi);
dP4deta =  0.25*(1-xi);
J(1,1) = v1(1)*dP1dxi ...
       + v2(1)*dP2dxi ...
       + v3(1)*dP3dxi ...
       + v4(1)*dP4dxi;

J(1,2) = v1(1)*dP1deta ...
       + v2(1)*dP2deta ...
       + v3(1)*dP3deta ...
       + v4(1)*dP4deta;

J(2,1) = v1(2)*dP1dxi ...
       + v2(2)*dP2dxi ...
       + v3(2)*dP3dxi ...
       + v4(2)*dP4dxi;

J(2,2) = v1(2)*dP1deta ...
       + v2(2)*dP2deta ...
       + v3(2)*dP3deta ...
       + v4(2)*dP4deta;
end


function J = quad_jacobian2(x,y,xi,eta)
J = zeros(2,2);
J(1,1) = ( (x(2)-x(1))*(1-eta) ...
         + (x(3)-x(4))*eta );
J(2,1) = ( (y(2)-y(1))*(1-eta) ...
         + (y(3)-y(4))*eta );
J(1,2) = ( (x(4)-x(1))*(1-xi) ...
         + (x(3)-x(2))*xi );
J(2,2) = ( (y(4)-y(1))*(1-xi) ...
         + (y(3)-y(2))*xi );
J = 0.5*J;
end

function c = quad_centroid(x,y)
% P1 = 0.25*(1-xi).*(1-eta);
% P2 = 0.25*(1+xi).*(1-eta);
% P3 = 0.25*(1+xi).*(1+eta);
% P4 = 0.25*(1-xi).*(1+eta);
c(1,1) = 0.25*( x(1) + x(2) + x(3) + x(4) );
c(2,1) = 0.25*( y(1) + y(2) + y(3) + y(4) );
end

function c = quad_centroid_2(x,y,npts)
X1 = [x(1),x(2);x(4),x(3)];
X2 = [y(1),y(2);y(4),y(3)];
X3 = [0,0;0,0];
quad_ref = quad_t.create_quad_ref_2D(npts);
quad_phys = quad_t.map_quad_ref_2D(X1,X2,X3,quad_ref,lagrange_interpolant());
volume = quad_phys.integrate(ones(1,quad_phys.n_quad));
c = quad_phys.integrate(quad_phys.quad_pts(1:2,:))/volume;
end

function [zeta,theta,phi,rho] = parameters_from_jacobian(A)
L1 = vmag( A(:,1) );
L2 = vmag( A(:,2) );
theta  = atan2( A(2,1), A(1,1) );
theta_ = atan2( A(2,2), A(1,2) );
phi    = theta_ - theta;
rho    = L2/L1;
zeta   = L1*L2;
end

function W = parameterized_2D_matrix(zeta,theta,phi,rho)
W = zeros(2,2);

W(1,1) = (1/sqrt(rho))*cos(theta);
W(2,1) = (1/sqrt(rho))*sin(theta);
W(1,2) =    sqrt(rho) *cos(theta+phi);
W(2,2) =    sqrt(rho) *sin(theta+phi);
W = sqrt(zeta)*W;
end

function [rt_zeta,R,Q,D,S,U] = jacobian_decomp(A)
[zeta,theta,phi,rho] = parameters_from_jacobian(A);
rt_zeta = sqrt(zeta);
rt_rho = sqrt(rho);
R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
Q = [1, cos(phi); 0, sin(phi)];
D = [1/rt_rho,0;0,rt_rho];
S = Q*D;
U = rt_zeta*S;
end