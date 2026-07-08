%% Grid Metric Prototyping v5 (07/07/2026)
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


%% make sure random quad is convex, and sort the points to be counterclockwise
xAR2 = -1;
sz = 5;
while xAR2<=0
    [x_vec,y_vec] = sort_points_counterclockwise(sz*rand([4,1]),sz*rand([4,1]));
    xAR2 = aspect_ratio_2(x_vec,y_vec);
end

% x_vec = [0,1,1,0];
% y_vec = [0,0,1,1];

%% set a point in computational space to evaluate at 
xi = 0; eta = 0;
xi_vec  = [-1,1,1,-1];
eta_vec = [-1,-1,1,1];


%% compute centroid
c1 = polygon_centroid(x_vec,y_vec);
c2 = [mean(x_vec);mean(y_vec)];

% find where the area centroid is in parametric space
xe_c1 = inverse_bilinear_map([x_vec(1),y_vec(1)],[x_vec(2),y_vec(2)],...
                            [x_vec(3),y_vec(3)],[x_vec(4),y_vec(4)],c1);
xe_c1 = xe_c1*2 - 1; % map to [-1,1]

hold on
axis equal
% plot(xi,eta,'ko')
% plot(xe_c1(1),xe_c1(2),'rx')
% [zeta,theta,phi,rho] = parameters_from_jacobian(A);
% W = parameterized_2D_matrix(zeta,theta,phi,rho);
% [rt_zeta,R,Q,D,S,U] = jacobian_decomp(A);

%% compute Jacobian at specified [xi,eta]
A = quad_jacobian2(x_vec,y_vec,xi,eta);

%% plot the cell in physical space, along with the basis (tangent) vectors
plot([x_vec(:);x_vec(1)],[y_vec(:);y_vec(1)])
[xp,yp] = quad_interp(x_vec,y_vec,xi,eta);
plot([xp,xp+A(1,1)],[yp,yp+A(2,1)],'r')
plot([xp,xp+A(1,2)],[yp,yp+A(2,2)],'b')

%% translate and apply the inverse mapping to all points
tmp = A\[(x_vec-c1(1)).';(y_vec-c1(2)).'];
x2 = tmp(1,:).';
y2 = tmp(2,:).';

tmp = A\[(x_vec-c2(1)).';(y_vec-c2(2)).'];
x2_ = tmp(1,:).';
y2_ = tmp(2,:).';

%% plot the reference square
plot([xi_vec(:);xi_vec(1)],[eta_vec(:);eta_vec(1)],'k')

%% plot the transformed points (should be closer to the reference square)
plot([x2;x2(1)],[y2;y2(1)],'--')
plot([x2_;x2_(1)],[y2_;y2_(1)],'.-')

%% plot the basis (tangent) vectors in the transformed space
A2 = quad_jacobian2(x2,y2,xi,eta);
[xp,yp] = quad_interp(x2,y2,xi,eta);
plot([xp,xp+A2(1,1)],[yp,yp+A2(2,1)],'r--')
plot([xp,xp+A2(1,2)],[yp,yp+A2(2,2)],'b--')

%% Now, find the optimal (L2) affine transform
%% ( approximate the transformation [xi_vec,eta_vec]=>[x_vec,y_vec] as affine )
% https://www3.cs.stonybrook.edu/~mueller/teaching/cse616/p17-rezk-salama.pdf
Aff = approximate_transform_as_affine(xi_vec,eta_vec,x_vec,y_vec);
%% apply the inverse transform and plot it
[x3,y3] = apply_inverse_affine_transform(Aff,x_vec,y_vec);

% ... it seems to lie on the inverse jacobian transform when evaluated at the
%     vertex centroid
plot([x3(:);x3(1)],[y3(:);y3(1)],'m--')

function A = approximate_transform_as_affine(x1,y1,x2,y2)
M = zeros(3);
for i = 1:4
    vec1 = [x1(i),y1(i),1];
    M = M + vec1(:)*(vec1(:))';
end
MAT = zeros(3);
for i = 1:4
    vec1 = [x1(i),y1(i),1];
    vec2 = [x2(i),y2(i),1];
    MAT = MAT + vec1(:)*(vec2(:))';
end
A = (M\MAT)';
end


function [x2,y2] = apply_affine_transform(A,x1,y1)
tmp = A*[x1(:).';y1(:).';1+0*x1(:).'];
x2 = tmp(1,:).';
y2 = tmp(2,:).';
end

function [x2,y2] = apply_inverse_affine_transform(A,x1,y1)
tmp = A\[x1(:).';y1(:).';1+0*x1(:).'];
x2 = tmp(1,:).';
y2 = tmp(2,:).';
end



function val = vmag(v)
val = sqrt( sum(v.^2) );
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

function [xp,yp] = quad_interp(x,y,xi,eta)
P(1,1) = (1-xi).*(1-eta);
P(2,1) = (1+xi).*(1-eta);
P(3,1) = (1+xi).*(1+eta);
P(4,1) = (1-xi).*(1+eta);
xp = 0.25*dot(P,x);
yp = 0.25*dot(P,y);
end

function J = quad_jacobian2(x,y,xi,eta)
J = zeros(2,2);
J(1,1) =  (x(2)-x(1))*(1-eta) ...
        + (x(3)-x(4))*(1+eta);
J(2,1) =  (y(2)-y(1))*(1-eta) ...
        + (y(3)-y(4))*(1+eta);
J(1,2) =  (x(4)-x(1))*(1-xi) ...
        + (x(3)-x(2))*(1+xi);
J(2,2) =  (y(4)-y(1))*(1-xi) ...
        + (y(3)-y(2))*(1+xi);
J = 0.25*J;
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


function out=argEval(n,fh, varargin)
%Get the n-th output of fh(varargin{:})
   [argsout{1:n}]=fh(varargin{:});
   
   out=argsout{n};
end

function out=scalarEval(fh,idx,varargin)
out1 = fh(varargin{:});
out = out1(idx);
end

function out=vectorInEval(fh,vec,varargin)
out = fh(vec{:},varargin{:});
end