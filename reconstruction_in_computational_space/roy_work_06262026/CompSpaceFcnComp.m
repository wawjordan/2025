clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up computational space and transformation matrices %
scale = 1.0; % NOTE: initial testing showed a scale near unity was best
xsi = [-1.;0.;1.;-1.;0.;1.;-1.;0.;1.];
eta = [-1.;-1.;-1.;0.;0.;0.;1.;1.;1.];
xsi = xsi*scale;  
eta = eta*scale;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set up grid in physical space (NOMINAL VALUES)
% ds = 1.0;        % initial streamwise spacing
% Ss = 1.0;        % streamwise stretching rate
% dn = 1.0;        % initial normal spacing
% Sn = 1.0;        % normal stretching rate
% theta = 0;       % grid angle in degrees
% xoffset = 1;     % x offset
% yoffset = 1;     % y offset
% betadeg = -30.0; % rotation angle in degrees
% beta = (betadeg./180.)*pi;

% Set up grid in physical space
ds = 1.0;        % initial streamwise spacing
Ss = 1.4;        % streamwise stretching rate
dn = 0.1;        % initial normal spacing
Sn = 1.4;        % normal stretching rate
theta = 25;      % grid angle in degrees
xoffset = 1;     % x offset
yoffset = 1;     % y offset
betadeg = -30.0; % rotation angle in degrees
beta = (betadeg./180.)*pi;


% Set up range of geometric parameters
dsvec = exp(linspace(log(0.01),log(100),9));
Ssvec = exp(linspace(log(0.01),log(100),9));
dnvec = exp(linspace(log(0.01),log(100),9));
Snvec = exp(linspace(log(0.01),log(100),9));
thetavec = linspace(0,320,9);
betadegvec = linspace(0,320,9);

% Loop over geometric parameters
for np = 1:9  % Only need to loop if sweeping in a parameter

    % NOTE: uncomment one of these below to sweep in that parameter
%     ds = dsvec(np);
%     Ss = Ssvec(np);
%     dn = dnvec(np);
%     Sn = Snvec(np);
%     theta = thetavec(np);
%     betadeg = betadegvec(np);
%     beta = (betadeg./180.)*pi;
    
    % Reset 9x9 computational transformation matrix
    for k = 1:9
        M(k,:) = [1,xsi(k),eta(k),xsi(k)*eta(k),xsi(k)^2,eta(k)^2,xsi(k)^2*eta(k)^2,xsi(k)^2*eta(k),xsi(k)*eta(k)^2];
    end
    M;
    condM = cond(M);
    Minv = inv(M);
    cond(Minv);
    
    % Set up grid in physical space
    x(1,1) = 0.;
    x(2,1) = x(1) + ds;
    x(3,1) = x(2) + Ss*ds*cosd(theta);
    x(4,1) = x(1);
    x(5,1) = x(2) + dn*tand(theta/2);
    x(6,1) = x(5) + Ss*(x(5)-x(4))*cosd(theta);
    x(7,1) = x(1);
    x(8,1) = x(2) + (dn+dn*Sn)*tand(theta/2);
    x(9,1) = x(8) + Ss*(x(8)-x(7))*cosd(theta);
    
    %x = transpose(x);
    
    y(1,1) = 0.;
    y(2,1) = y(1);
    y(3,1) = y(2) - Ss*ds*sind(theta);
    y(4,1) = y(1) + dn;
    y(5,1) = y(4);
    y(6,1) = y(5) - Ss*ds*sind(theta);
    y(7,1) = y(4) + Sn*dn;
    y(8,1) = y(7);
    y(9,1) = y(8) - Ss*ds*sind(theta);
    
    % Make y a column vector
    %y = transpose(y);
    
    % Offsets in x and y
    x = x + xoffset;
    y = y + yoffset;
    
    % Rotate by angle beta
    r = sqrt(x.^2 + y.^2);
    gamma = atan(y./x);
    delta = gamma + beta;
    y = r.*sin(delta);
    x = r.*cos(delta);
    
    % Reshape into matrices for plotting purposes
    xx = reshape(x,[3,3]);
    yy = reshape(y,[3,3]);
    zz = [1,1,1;1,1,1;1,1,1];
    
    % Plot the 2D grid
    figure(1);
    surf(xx, yy, zz);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Grid');
    daspect([1 1 1]);
    view(2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for the coefficients of a vector = [a1, a2, ..., a9]^T
    
    format long
    
    % for the x values
    axvec = Minv*x;
    
    % for the y values
    ayvec = Minv*y;
    
    % Checks (should return zero within machine precision)
    checkx = M*axvec - x;
    checky = M*ayvec - y;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define a function to be reconstructed in computational space
    funct = functcomp(xsi,eta);
    
    % Make into a matrix for plotting
    functmat = reshape(funct,[3,3]);
    
    % plot the function
    figure(2);
    surf(xx, yy, functmat);
    xlabel('X');
    ylabel('Y');
    zlabel('funct');
    title('Function');
    daspect([1 1 1]);
    %view(2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reconstruct in computational space
    
    bvec = Minv*funct;
    
    % Check that function is reconstructed correctly at nodes (w/in machine
    % precision)
    
    err = funct - M*bvec;
    
    
    % Reshape xsi and eta vectors into matrices
    xsimat = reshape(xsi,[3,3]);
    etamat = reshape(eta,[3,3]);
    
    % How well is function reconstructed at cell centers in comp space?
    
    % Calculate computational and physical cell centers in matrix form
    for i = 1:2
        for j = 1:2
            % Computational space
            xsicc(i,j) = 0.25 * (xsimat(i,j) + xsimat(i+1,j) + xsimat(i,j+1) + xsimat(i+1,j+1) );
            etacc(i,j) = 0.25 * (etamat(i,j) + etamat(i+1,j) + etamat(i,j+1) + etamat(i+1,j+1) );
            % Physical space
            xxcc(i,j) = 0.25 * (xx(i,j) + xx(i+1,j) + xx(i,j+1) + xx(i+1,j+1) );
            yycc(i,j) = 0.25 * (yy(i,j) + yy(i+1,j) + yy(i,j+1) + yy(i+1,j+1) );
        end
    end
    
    xsicc;
    etacc;
    xxcc;
    yycc;
    
    % Convert cell centers to vector form
    xsiccvec = xsicc(:);
    etaccvec = etacc(:);
    xcc = xxcc(:);
    ycc = yycc(:);
    
    % Reconstruct at computational cell centers
    for k = 1:4
        xsi2 = xsiccvec(k);
        eta2 = etaccvec(k);
        Mrecon(k,:) = [1,xsi2,eta2,xsi2*eta2,xsi2^2,eta2^2,xsi2^2*eta2^2,xsi2^2*eta2,xsi2*eta2^2];
    end
    
    % Do cell centers in computational space equal cell centers in physical
    % space? ==> should only be so for uniform meshes (maybe)
    
    xtran = Mrecon*axvec;
    ytran = Mrecon*ayvec;
    
    errtranx = xtran - xcc;
    errtrany = ytran - ycc;
    
    errtranxnorm = norm(errtranx);
    errtranynorm = norm(errtrany);
    
    Mrecon;
    
    % Error in reconstruction at computational grid cell centers
    frecon = Mrecon*bvec;
    functcc = functcomp(xsiccvec,etaccvec);
    
    err = frecon - functcc;
    errcompccrecon = norm(err);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find physical cell center locations in computational space (iteratively)
    
    nmax = 1000; % max number of iterations
    for k = 1:4
        % Target x and y are physical cell centers
        xtarg = xcc(k);
        ytarg = ycc(k);
        % Start guess at computational cell center
        xsin = xsiccvec(k);
        etan = etaccvec(k);
        for n = 1:nmax
            Mreconn = [1,xsin,etan,xsin*etan,xsin^2,etan^2,xsin^2*etan^2,xsin^2*etan,xsin*etan^2];
            Rx = Mreconn*axvec - xtarg;
            Ry = Mreconn*ayvec - ytarg;
            %disp([num2str(n),',  ',num2str(Rx),',  ',num2str(Ry)]);
            dRxdxsi = axvec(2) + axvec(4)*etan + 2.*axvec(5)*xsin + 2.*axvec(7)*xsin*etan^2 + 2.*axvec(8)*xsin*etan + axvec(9)*etan^2;
            dRxdeta = axvec(3) + axvec(4)*xsin + 2.*axvec(6)*etan + 2.*axvec(7)*xsin^2*etan + axvec(8)*xsin^2 + 2.*axvec(9)*xsin*etan;
            dRydxsi = ayvec(2) + ayvec(4)*etan + 2.*ayvec(5)*xsin + 2.*ayvec(7)*xsin*etan^2 + 2.*ayvec(8)*xsin*etan + ayvec(9)*etan^2;
            dRydeta = ayvec(3) + ayvec(4)*xsin + 2.*ayvec(6)*etan + 2.*ayvec(7)*xsin^2*etan + ayvec(8)*xsin^2 + 2.*ayvec(9)*xsin*etan;
            detJ = dRxdxsi*dRydeta - dRxdeta*dRydxsi;
            xsin = xsin - (Rx*dRydeta - Ry*dRxdeta)/detJ;
            etan = etan - (Ry*dRxdxsi - Rx*dRydxsi)/detJ;
        end
        if (Rx > 1E-10) || (Ry > 1E-10)
            disp([num2str(k),',  ',num2str(Rx),',  ',num2str(Ry)]);
            error('Error, iterations did not converge. Try increasing nmax.')
        end
        xsitarg(k,1) = xsin;
        etatarg(k,1) = etan;
    end
    
    % Make target xsi and eta vectors column vectors
    %xsitarg = transpose(xsitarg);
    %etatarg = transpose(etatarg);
    
    
    % Calculate error in computational space reconstruction at physical cell
    % centers (found via iterative scheme)
    
    % Reconstruct at cell physical grid centers
    for k = 1:4
        xsi2 = xsitarg(k);
        eta2 = etatarg(k);
        Mrecon(k,:) = [1,xsi2,eta2,xsi2*eta2,xsi2^2,eta2^2,xsi2^2*eta2^2,xsi2^2*eta2,xsi2*eta2^2];
    end
    
    % Error in reconstruction using physical grid cell centers
    frecon = Mrecon*bvec;
    functcc = functcomp(xsitarg,etatarg);
    err = frecon - functcc;
    errphysccrecon = norm(err);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOW, need to try reconstructing in physical space only
    
    for k = 1:9
        G(k,:) = [1,x(k),y(k),x(k)*y(k),x(k)^2,y(k)^2,x(k)^2*y(k)^2,x(k)^2*y(k),x(k)*y(k)^2];
    end
    
    G;
    condG = cond(G);
    Ginv = inv(G);
    cond(Ginv);
    
    cvec = Ginv*funct;
    
    funct - G*cvec;  % check to make sure nodal functions is reconstructed correctly w/in machine precision
    
    
    % Calculate error in physical space reconstruction at physical cell centers
    % Reconstruct at cell centers
    for k = 1:4
        Grecon(k,:) = [1,xcc(k),ycc(k),xcc(k)*ycc(k),xcc(k)^2,ycc(k)^2,xcc(k)^2*ycc(k)^2,xcc(k)^2*ycc(k),xcc(k)*ycc(k)^2];
    end
    
    % Error in reconstruction using physical grid cell centers
    freconphys = Grecon*cvec;
    err = freconphys - functcc;
    errphysccreconphys = norm(err);
    
    format long
disp([num2str(ds),'  ',num2str(Ss),'  ',num2str(dn),'  ',num2str(Sn) ...
        ,'  ',num2str(theta),'  ',num2str(betadeg),'  ' ...
        ,num2str(errcompccrecon),'  ',num2str(errphysccrecon) ...
        ,'  ',num2str(errphysccreconphys)...
        ,'  ',num2str(condM),'  ',num2str(condG)])

end


%************************************************************************
function f = functcomp(xsi,eta)  
f = 0.5 + 2*sin( (pi *(1.*xsi.*eta + 0.7.*xsi + 0.3.*eta + 0.25) ) / 20); 

end

