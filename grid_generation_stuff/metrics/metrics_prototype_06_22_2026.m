%% Grid Metric Prototyping
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

v1 = [0,0];
v2 = [1,0];
v3 = [1,2];
v4 = [0,1];
[x,y] = sort_points_counterclockwise(rand([4,1]),rand([4,1]));
v1 = [x(1),y(1)];
v2 = [x(2),y(2)];
v3 = [x(3),y(3)];
v4 = [x(4),y(4)];

s1  = skewness_1(v1,v2,v3,v4);
s1xy = skewness_1_xy(x,y);
s2  = skewness_2(v1,v2,v3,v4);
s2xy = skewness_2xy(x,y);
T   = taper(v1,v2,v3,v4);
T_xy= taper_xy(x,y);
AR1 = aspect_ratio_1(v1,v2,v3,v4);
AR1xy = aspect_ratio_1_xy(x,y);
AR2 = aspect_ratio_2(v1,v2,v3,v4);
AR2xy = aspect_ratio_2xy(x,y);
AR3 = aspect_ratio_3(v1,v2,v3,v4);
AR3xy = aspect_ratio_3xy(x,y);

% A1 = quad_jacobian2(v1,v2,v3,v4,0,0)
j1 = jac_ind_1(v1,v2,v3,v4);
j2 = jac_ind_1_knupp(v1,v2,v3,v4);

j2xy = jac_ind_1_knupp_xy(x,y);

xy = [v1;v2;v3;v4;v1];
plot(xy(:,1),xy(:,2))
axis equal


function val = vmag(v)
val = sqrt( sum(v.^2) );
end

function [x1,x2] = quad_principal_axes(v1,v2,v3,v4)
% v4 +-<-+ v3
%    |   |
% v1 +->-+ v2
x1 = ( v2 - v1 ) + (v3 - v4 );
x2 = ( v3 - v2 ) + (v4 - v1 );
end

function val = skewness_1(v1,v2,v3,v4)
[x1,x2] = quad_principal_axes(v1,v2,v3,v4);
val = abs( dot( x1./vmag(x1), x2./vmag(x2) ) );
end

function val = skewness_1_xy(x,y)
x1   = ( x(2) - x(1) ) + (x(3) - x(4) );
y1   = ( y(2) - y(1) ) + (y(3) - y(4) );
x2   = ( x(3) - x(2) ) + (x(4) - x(1) );
y2   = ( y(3) - y(2) ) + (y(4) - y(1) );
m1m2 = sqrt( x1^2 + y1^2 ) * sqrt( x2^2 + y2^2 );
val  = abs( (x1*x2+y1*y2)/m1m2 );
end

function val = skewness_2(v1,v2,v3,v4)
v41 = v4 - v1;
v21 = v2 - v1;
v32 = v3 - v2;
v43 = v4 - v3;
xs1 = 1/sqrt( 1 - ( dot( v41, v21 )/(vmag(v41)*vmag(v21)) )^2 );
xs2 = 1/sqrt( 1 - ( dot( v21, v32 )/(vmag(v21)*vmag(v32)) )^2 );
xs3 = 1/sqrt( 1 - ( dot( v32, v43 )/(vmag(v32)*vmag(v43)) )^2 );
xs4 = 1/sqrt( 1 - ( dot( v43, v41 )/(vmag(v43)*vmag(v41)) )^2 );
val = 4/(xs1 +xs2+xs3+xs4);
end

function val = skewness_2xy(x,y)
% v41 = v4 - v1;
% v21 = v2 - v1;
% v32 = v3 - v2;
% v43 = v4 - v3;
% xs1 = 1/sqrt( 1 - ( dot( v41, v21 )/(vmag(v41)*vmag(v21)) )^2 );
% xs2 = 1/sqrt( 1 - ( dot( v21, v32 )/(vmag(v21)*vmag(v32)) )^2 );
% xs3 = 1/sqrt( 1 - ( dot( v32, v43 )/(vmag(v32)*vmag(v43)) )^2 );
% xs4 = 1/sqrt( 1 - ( dot( v43, v41 )/(vmag(v43)*vmag(v41)) )^2 );
% val = 4/(xs1 +xs2+xs3+xs4);

v0x = x(1)-x(4);
v0y = y(1)-y(4);
m0  = sqrt( v0x^2 + v0y^2 );

v2x = v0x;
v2y = v0y;
m2  = m0;

val = 0;
for i = 2:4
    v1x = x(i)-x(i-1);
    v1y = y(i)-y(i-1);
    m1  = sqrt( v1x^2 + v1y^2 );
    val = val + 1/sqrt( 1 - ( (v1x*v2x + v1y*v2y)/(m1*m2) )^2 );
    m2  = m1;
    v2x = v1x;
    v2y = v1y;
end
val = val + 1/sqrt( 1 - ( (v1x*v0x + v1y*v0y)/(m1*m0) )^2 );
val = 4/val;
end

function val = taper(v1,v2,v3,v4)
[x1,x2] = quad_principal_axes(v1,v2,v3,v4);
x12 = ( v1 - v2 ) + (v3 - v4 );
val = vmag(x12)/min( vmag(x1), vmag(x2) );
end

function val = taper_xy(x,y)
m1  = sqrt( ( ( x(2) - x(1) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(2) - y(1) ) + (y(3) - y(4) ) )^2 );
m2  = sqrt( ( ( x(3) - x(2) ) + (x(4) - x(1) ) )^2 ...
          + ( ( y(3) - y(2) ) + (y(4) - y(1) ) )^2 );
m12 = sqrt( ( ( x(1) - x(2) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(1) - y(2) ) + (y(3) - y(4) ) )^2 );
val = m12/min( m1, m2 );
end

function val = aspect_ratio_1(v1,v2,v3,v4)
[x1,x2] = quad_principal_axes(v1,v2,v3,v4);
val = vmag(x1)/vmag(x2);
val = max( val, 1/val );
end

function val = aspect_ratio_1_xy(x,y)
m1  = sqrt( ( ( x(2) - x(1) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(2) - y(1) ) + (y(3) - y(4) ) )^2 );
m2  = sqrt( ( ( x(3) - x(2) ) + (x(4) - x(1) ) )^2 ...
          + ( ( y(3) - y(2) ) + (y(4) - y(1) ) )^2 );
val = m1/m2;
val = max( val, 1/val );
end

function val = aspect_ratio_2(v1,v2,v3,v4)
h1 = vmag(v2-v1);
h2 = vmag(v3-v2);
h3 = vmag(v4-v3);
h4 = vmag(v1-v4);
d1 = vmag(v3-v1);
d2 = vmag(v4-v2);
% [(x2-x1)(y3-y1)-(x3-x1)(y2-y1)]


A1 = 0.5*( ( v2(1)-v1(1) )*( v3(2)-v1(2) ) ...
         - ( v3(1)-v1(1) )*( v2(2)-v1(2) ) );
A2 = 0.5*( ( v3(1)-v2(1) )*( v4(2)-v2(2) ) ...
         - ( v4(1)-v2(1) )*( v3(2)-v2(2) ) );
A3 = 0.5*( ( v4(1)-v3(1) )*( v1(2)-v3(2) ) ...
         - ( v1(1)-v3(1) )*( v4(2)-v3(2) ) );
A4 = 0.5*( ( v2(1)-v1(1) )*( v4(2)-v1(2) ) ...
         - ( v4(1)-v1(1) )*( v2(2)-v1(2) ) );

% A1_1 = -polygon_area([v1(1),v2(1),v3(1)],[v1(2),v2(2),v3(2)]);
% A2_1 = -polygon_area([v2(1),v3(1),v4(1)],[v2(2),v3(2),v4(2)]);
% A3_1 = -polygon_area([v3(1),v4(1),v1(1)],[v3(2),v4(2),v1(2)]);
% A4_1 = -polygon_area([v4(1),v1(1),v2(1)],[v4(2),v1(2),v2(2)]);
hmax = max([h1,h2,h3,h4]);
alpha = sqrt(2)/8;
val = alpha*vmag([h1,h2,h3,h4])*max([hmax,d1,d2])/min([A1,A2,A3,A4]);
end


function val = aspect_ratio_2xy(x,y)
% minA =            (x(1)-x(4))*(y(2)-y(4))-(x(2)-x(4))*(y(1)-y(4));
% minA = min( minA, (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)) );
% minA = min( minA, (x(3)-x(2))*(y(4)-y(2))-(x(4)-x(2))*(y(3)-y(2)) );
% minA = min( minA, (x(4)-x(3))*(y(1)-y(3))-(x(1)-x(3))*(y(4)-y(3)) );
% minA = 0.5*minA;
% hmax = 0;
% hmag = 0;
% h2   = (x(1)-x(4))^2 + (y(1)-y(4))^2;
% hmag = hmag + h2;
% hmax = max( hmax, sqrt(h2) );
% 
% h2   = (x(2)-x(1))^2 + (y(2)-y(1))^2;
% hmag = hmag + h2;
% hmax = max( hmax, sqrt(h2) );
% 
% h2   = (x(3)-x(2))^2 + (y(3)-y(2))^2;
% hmag = hmag + h2;
% hmax = max( hmax, sqrt(h2) );
% 
% h2   = (x(4)-x(3))^2 + (y(4)-y(3))^2;
% hmag = hmag + h2;
% hmax = max( hmax, sqrt(h2) );
% 
% hmag = sqrt( hmag );
% 
% hmax = max( hmax, sqrt( (x(3)-x(1))^2 + (y(3)-y(1))^2 ) );
% hmax = max( hmax, sqrt( (x(4)-x(2))^2 + (y(4)-y(2))^2 ) );
% 
% val = (sqrt(2)/8)*hmag*hmax/minA;

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
val = 0.25*sqrt(2*hmag)*hmax/minAx2;
end

function val = aspect_ratio_3(v1,v2,v3,v4)
v41 = v4 - v1;
v21 = v2 - v1;
v32 = v3 - v2;
v43 = v4 - v3;
a1 = acos( ( dot( v41, v21 )/(vmag(v41)*vmag(v21)) ) );
a2 = acos( ( dot(-v21, v32 )/(vmag(v21)*vmag(v32)) ) );
a3 = acos( ( dot(-v32, v43 )/(vmag(v32)*vmag(v43)) ) );
a4 = acos( ( dot(-v43,-v41 )/(vmag(v43)*vmag(v41)) ) );
alpha1 = 0.5*(a1+a3);
alpha2 = 0.5*(a2+a4);
val = max( 1/(1-abs(1-(2*alpha1)/pi) ), 1/(1-abs(1-(2*alpha2)/pi) ) );
end

function val = aspect_ratio_3xy(x,y)

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
val = 1/(1-abs(1-(2*alpha)/pi) );
end

% function val = oddy(v1,v2,v3,v4)
% val = 0;
% h1 = vmag(v2-v1);
% h2 = vmag(v3-v2);
% h3 = vmag(v4-v3);
% h4 = vmag(v1-v4);
% end
% [coefs,d_order] = finite_difference_coefs_v1(stencil,d_order);
% function [x_xi,y_xi,x_eta,y_eta] = coord_diffs_C_grid(x,y,order)
% x_xi  = zeros(size(x));
% y_xi  = zeros(size(x));
% x_eta = zeros(size(x));
% y_eta = zeros(size(x));
% end

function val = jac_ind_1(v1,v2,v3,v4)
den = 0;
cnt = 0;
% for j = -1:2:1
%     for i = -1:2:1
for j = 0:1
    for i = 0:1
        cnt = cnt + 1;
        A = quad_jacobian2(v1,v2,v3,v4,i,j);
        A2 = quad_jacobian_knupp(v1,v2,v3,v4,cnt);
        A = A'*A;
        den = den + sum(eig(A))/det(A);
    end
end
val = 8/den;
end

function val = jac_ind_1_knupp(v1,v2,v3,v4)
den = 0;
for n = 1:4
    A = quad_jacobian_knupp(v1,v2,v3,v4,n);
    A = A'*A;
    den = den + sum(eig(A))/det(A);
end
val = 8/den;
end

function val = jac_ind_1_knupp_xy(x,y)
den = 0;
for n = 1:4
    A = quad_jacobian_knupp_xy(x,y,n);
    A = A'*A;
    den = den + sum(eig(A))/det(A);
end
val = 8/den;
end

function J = quad_jacobian(v1,v2,v3,v4,xi,eta)
J = zeros(2,2);

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

J(2,1) = v1(1)*dP1deta ...
       + v2(1)*dP2deta ...
       + v3(1)*dP3deta ...
       + v4(1)*dP4deta;

J(1,2) = v1(2)*dP1dxi ...
       + v2(2)*dP2dxi ...
       + v3(2)*dP3dxi ...
       + v4(2)*dP4dxi;

J(2,2) = v1(2)*dP1deta ...
       + v2(2)*dP2deta ...
       + v3(2)*dP3deta ...
       + v4(2)*dP4deta;
end


function J = quad_jacobian2(v1,v2,v3,v4,xi,eta)
J = zeros(2,2);
% J(1,:) = 0.25*( (v2-v1)*(1-eta) ...
%               + (v3-v4)*(1+eta) );
% J(2,:) = 0.25*( (v4-v1)*(1-xi) ...
%               + (v3-v2)*(1+xi) );
J(1,:) = ( (v2-v1)*(1-eta) ...
         + (v3-v4)*eta );
J(2,:) = ( (v4-v1)*(1-xi) ...
         + (v3-v2)*xi );
end

function J = quad_jacobian_knupp(v1,v2,v3,v4,k)
J = zeros(2,2);
xy = [v1(:).';v2(:).';v3(:).';v4(:).'];
J(1,1) = xy(mod(k-1+1,4)+1,1) - xy(mod(k-1,4)+1,1);
J(1,2) = xy(mod(k-1+3,4)+1,1) - xy(mod(k-1,4)+1,1);
J(2,1) = xy(mod(k-1+1,4)+1,2) - xy(mod(k-1,4)+1,2);
J(2,2) = xy(mod(k-1+3,4)+1,2) - xy(mod(k-1,4)+1,2);
end

function J = quad_jacobian_knupp_xy(x,y,k)
J = zeros(2,2);
kp1  = mod(k  ,4)+1;
kp3  = mod(k+2,4)+1;
J(1,1) = x(kp1) - x(k);
J(1,2) = x(kp3) - x(k);
J(2,1) = y(kp1) - y(k);
J(2,2) = y(kp3) - y(k);
end
