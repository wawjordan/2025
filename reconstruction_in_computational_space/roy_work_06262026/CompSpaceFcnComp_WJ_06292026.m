%% Edited Version of Script from Dr. Roy (06/29/2026)
clc; clear; close all;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parent_dir_str = '2025';
% path_parts = regexp(mfilename('fullpath'), filesep, 'split');
% path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
% parent_dir = fullfile(path_parts{1:path_idx});
% addpath(genpath(parent_dir));
% clear parent_dir_str path_idx path_parts
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;


% Set up grid in physical space
in = struct();
in.ds       =   1.0; % initial streamwise spacing
in.Ss       =   1.4; % streamwise stretching rate
in.dn       =   0.1; % initial normal spacing
in.Sn       =   1.4; % normal stretching rate
in.thetadeg =  25.0; % grid angle in degrees
in.xoffset  =   1.0; % x offset
in.yoffset  =   1.0; % y offset
in.betadeg  = -30.0; % rotation angle in degrees

in.degree   = 2;
in.scale    = 1.0;

xsi_vec = -1:1;
eta_vec = -1:1;

[XSI,ETA] = ndgrid(xsi_vec,eta_vec);

[x,y] = generate_grid( in.ds, in.Ss, in.dn, in.Sn, in.thetadeg, ...
                       in.xoffset, in.yoffset, in.betadeg, xsi_vec, eta_vec );

fvals = functcomp(XSI,ETA);

% Plot the 2D grid
figure(1);
surf(x, y, 0*x);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Grid');
daspect([1 1 1]);
view(2)

% plot the function
figure(2);
surf(x, y, fvals);
xlabel('X');
ylabel('Y');
zlabel('funct');
title('Function');
daspect([1 1 1]);

[M,exponents] = computational_transform_matrix({xsi_vec,eta_vec},in.degree,in.scale);

axvec = x(:)\M;
ayvec = y(:)\M;

bvec  = fvals(:)\M;

% get cell centers
xsicc = get_nodal_centroid(XSI);
etacc = get_nodal_centroid(ETA);
xcc = get_nodal_centroid(x);
ycc = get_nodal_centroid(y);

% Reconstruct at computational cell centers
for k = 1:4
    xsi2 = xsicc(k);
    eta2 = etacc(k);
    Mrecon(k,:) = [1,xsi2,eta2,xsi2*eta2,xsi2^2,eta2^2,xsi2^2*eta2^2,xsi2^2*eta2,xsi2*eta2^2];
end

Mrecon2 = get_reconstructed_cell_center_matrix(exponents,{xsicc,etacc});


figure
function vals = get_reconstructed_cell_center_matrix(exponents,X)
Nc = numel(X{1});
n_terms = size(exponents,2);
vals = zeros(Nc,n_terms);

for n = 1:n_terms
    coef_vec = zeros(n_terms,1);
    coef_vec(n) = 1;
    tmp = reconstruct_point_vals(coef_vec,exponents,X);
    vals(:,n) = tmp(:);
end
end

function vals = reconstruct_point_vals(coefs,exponents,X)
n_dim = numel(X);
sz = size(X{1});
vals = zeros(sz);
nvals = prod(sz);

xvec = zeros(1,n_dim);
for i = 1:nvals
    for d = 1:n_dim
        xvec(d) = X{d}(i);
    end
    vals(i) = reconstruct_point_val(coefs,exponents,xvec);
end
end

function val = reconstruct_point_val(coefs,exponents,x)
val = 0;
n_dim = numel(x);
n_terms = numel(coefs);
for i = 1:n_terms
    val = val + coefs(i)*evaluate_monomial( n_dim, exponents, i, x );
end
end

function [Acc] = get_nodal_centroid(A)
[Ni,Nj] = size(A);
Acc = zeros(Ni-1,Nj-1);
for j =1:Nj-1
    for i = Ni-1
        Acc(i,j) = 0.25*( A(i,j) + A(i+1,j) + A(i,j+1) + A(i+1,j+1) );
    end
end
end


function [x,y] = generate_grid(ds,Ss,dn,Sn,thetadeg,xoffset,yoffset,betadeg,~,~)
% ds       - initial streamwise spacing
% Ss       - streamwise stretching rate
% dn       - initial normal spacing
% Sn       - normal stretching rate
% thetadeg - grid angle in degrees
% xoffset  - x offset
% yoffset  - y offset
% betadeg  - rotation angle in degrees
% xsi      - xi coordinate
% eta      - eta coordinate
theta = (thetadeg/180)*pi;
beta  = (betadeg/180)*pi;

% Nxsi = numel(xsi);
% Neta = numel(eta);
% 
% xsi0  = xsi(1);  xsiN = xsi(Nxsi); xsi1 = xsi(1); xsi2 = xsi(2);
% eta0  = eta(1);  etaN = xsi(Neta);
% 
% thetafun = @(xsi) (xsi>xsi2).*theta.*((xsi-xsi1)/(xsiN-xsi1));
% xfun1 = @(xsi) xoffset + cos(thetafun(xsi)).*ds*((Ss^(Nxsi-1)).^((xsi-xsi0)/(xsiN-xsi0)) - 1)/(Ss-1);
% yfun1 = @(eta) yoffset + dn*((Sn^(Neta-1)).^((eta-eta0)/(etaN-eta0)) - 1)/(Sn-1);
% 
% % generate the initial spacings
% x0 = xfun1(xsi);
% y0 = yfun1(eta);

Nxsi = 3;
Neta = 3;

x = zeros(Nxsi,Neta);
y = zeros(Nxsi,Neta);

x(1) = 0.;
x(2) = x(1) + ds;
x(3) = x(2) + Ss*ds*cos(theta);
x(4) = x(1);
x(5) = x(2) + dn*tan(theta/2);
x(6) = x(5) + Ss*(x(5)-x(4))*cos(theta);
x(7) = x(1);
x(8) = x(2) + (dn+dn*Sn)*tan(theta/2);
x(9) = x(8) + Ss*(x(8)-x(7))*cos(theta);

y(1) = 0.;
y(2) = y(1);
y(3) = y(2) - Ss*ds*sin(theta);
y(4) = y(1) + dn;
y(5) = y(4);
y(6) = y(5) - Ss*ds*sin(theta);
y(7) = y(4) + Sn*dn;
y(8) = y(7);
y(9) = y(8) - Ss*ds*sin(theta);

x = x + xoffset;
y = y + yoffset;
    
r = sqrt(x.^2 + y.^2);
gamma = atan(y./x);
delta = gamma + beta;
y = r.*sin(delta);
x = r.*cos(delta);

end





function [M,exponents] = computational_transform_matrix(X,degree,scale)
flag = false;
if ~isa(X,"cell")
    flag = true;
end
n_dim = numel(X);
if (n_dim>3)
    error('only dimensions 1-3 are supported')
end
for i = 1:n_dim
    if ~isa(X{i},'numeric')
        flag = true;
        break;
    elseif ~isvector(X{i})
        flag = true;
        break;
    end
end
if (flag)
    error('X must be a cell array of 1D monotonic vectors')
end
switch n_dim
    case(1)
        xvec = ( X{1}(:) ).';
    case(2)
        [x1,x2] = ndgrid( X{1}, X{2} );
        xvec = [x1(:),x2(:)].';
    case(3)
        [x1,x2,x3] = ndgrid( X{1}, X{2}, X{3} );
        xvec = [x1(:),x2(:),x3(:)].';
end
xvec = xvec*scale;
% [exponents,~,~] = get_exponents( n_dim, degree );
exponents = get_exponents_tp( n_dim, degree );
n_terms = size(exponents,2);
n_nodes = prod( cellfun(@numel,X) );

M = ones(n_nodes,n_terms);
for j = 1:n_terms
    for i = 1:n_nodes
        M(i,j) = evaluate_monomial( n_dim, exponents,j,xvec(:,i));
    end
end
end

function [val,coef] = evaluate_monomial( n_dim, exponents, term, x )
val = 1.0;
coef = 1;
for d = 1:n_dim
    for i = exponents(d,term):-1:1
        val = val*x(d);
        coef = coef * i;
    end
end
end

function [dval,dcoef,coef] = evaluate_monomial_derivative( n_dim,     ...
                                                           exponents, ...
                                                           term,      ...
                                                           x,         ...
                                                           order )
dval  = 0.0;
dcoef = 1;
coef  = 1;
if ( any( exponents(:,term) - order(:) < 0 ) ); return; end
dval = 1.0;
for d = 1:n_dim
    for i = exponents(d,term):-1:exponents(d,term)-order(d)+1
        dcoef = dcoef * i;
    end
    for i = exponents(d,term)-order(d):-1:1
        coef = coef * i;
        dval = dval*x(d);
    end
end
end

function exponents = get_exponents_tp( n_dim, degree )
n_terms   = (degree+1)^n_dim;
exponents = zeros(n_dim,n_terms);
nsub = (degree+1)*ones(1,n_dim);
for j = 1:n_terms
    exponents(:,j) = global2local(j,nsub)-1;
end
if (n_dim== 2 && degree==2)
    exponents = exponents(:,[1,2,4,5,3,7,9,6,8]);
end
% % 1: [0,0]  1
% % 2: [1,0]  2
% % 3: [0,1]  4
% % 4: [1,1]  5
% % 5: [2,0]  3
% % 6: [0,2]  7
% % 7: [2,2]  9
% % 8: [2,1]  6
% % 9: [1,2]  8
end

% function exponents = get_exponents_tp( n_dim, degree )
% max_degree   = degree*n_dim;
% [exponents,~,~] = get_exponents( n_dim, max_degree );
% exponents(:,any(exponents>degree,1)) = [];
% % 1: [0,0]  1
% % 2: [1,0]  2
% % 3: [0,1]  3
% % 4: [1,1]  5
% % 5: [2,0]  4
% % 6: [0,2]  6
% % 7: [2,2]  9
% % 8: [2,1]  7
% % 9: [1,2]  8
% end

function [exponents,idx,diff_idx] = get_exponents( n_dim, degree )
n_terms   = nchoosek( n_dim + degree, degree );
exponents = zeros(n_dim,n_terms);
idx       = zeros(degree+1,1);
diff_idx  = zeros(n_dim,n_terms);
cnt = 0;
for curr_total_degree = 0:degree
    nsub(1:n_dim) = curr_total_degree + 1;
    N_full_terms = (curr_total_degree+1)^n_dim;
    for j = 0:N_full_terms
        tmp_exp = global2local(j+1,nsub)-1;
        if ( sum(tmp_exp) == curr_total_degree )
            cnt = cnt + 1;
            exponents(:,cnt) = tmp_exp;
        end
    end
    idx(curr_total_degree+1) = cnt;
end
diff_idx(1:n_dim,1:n_terms) = -1;
if ( degree == 0 ); return; end
for j = 1:idx(degree)
    tmp_exp = exponents(:,j);
    curr_total_degree = sum(tmp_exp);
    cnt = 0;
    for i = idx(curr_total_degree+1)+1:idx(curr_total_degree+2)
        if ( sum( abs(exponents(:,i) - tmp_exp) ) == 1 )
            cnt = cnt +1;
            diff_idx(cnt,j) = i;
        end
    end
end
end

function iSub = global2local(iG,nSub)
nDims = numel(nSub);
iSub = zeros(1,nDims);
if (nDims==1)
    iSub(1) = iG;
    return
end
p = prod(nSub);
iGtmp = iG;
for i = nDims:-1:1
    p = fix( p/nSub(i) );
    iTmp = mod(iGtmp-1,p)+1;
    iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
    iGtmp = iTmp;
end
end

% function iG = local2global(iSub,nSub)
% iSub = iSub(:);
% nSub = nSub(:);
% nDims = numel(iSub);
% p = 1;
% iG = 1;
% for i = 1:nDims
%     iG = iG + ( iSub(i) - 1 )*p;
%     p = p*nSub(i);
% end
% end

function f = functcomp(xsi,eta)  
f = 0.5 + 2.0*sin( (pi *(1.0*xsi.*eta + 0.7*xsi + 0.3*eta + 0.25) ) / 20.0);
end