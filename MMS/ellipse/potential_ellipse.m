classdef potential_ellipse
    properties
        a(1,1)           double = 1.0
        b(1,1)           double = 1.0
        l(1,1)           double = 0.0
        R(1,1)           double = 1.0
        alpha(1,1)       double = 0.0
        circulation(1,1) double = 0.0
        rho_inf(1,1)     double = 1.0
        v_inf(1,1)       double = 68.0
        p_inf(1,1)       double = 100000.0
        gamma(1,1)       double = 1.4
        rho_ref(1,1)     double = 1.0
        a_ref(1,1)       double = 340.0
        l_ref(1,1)       double = 1.0
        l_ref_grid(1,1)  double = 1.0
    end
    methods
%% Constructor
        function this = potential_ellipse(varargin)
            validScalarNum        = @(x) isnumeric(x) && isscalar(x);
            validScalarNonNegNum  = @(x) validScalarNum(x) && (x >= 0);
            validScalarPosNum     = @(x) validScalarNum(x) && (x > 0);
            p = inputParser;
            p.addOptional('a',1.0,validScalarNonNegNum);
            p.addOptional('b',1.0,validScalarNonNegNum);
            p.addOptional('alpha',0.0,validScalarNum)
            p.addOptional('circulation',0.0,validScalarNum)
            p.addOptional('rho_inf',1.0,validScalarPosNum)
            p.addOptional('v_inf',68.0,validScalarPosNum)
            p.addOptional('p_inf',100000.0,validScalarPosNum)
            p.addOptional('gamma',1.4,validScalarPosNum)
            p.addOptional('rho_ref',1.0,validScalarPosNum)
            p.addOptional('a_ref',340.0,validScalarPosNum)
            p.addOptional('l_ref',1.0,validScalarPosNum)
            p.addOptional('l_ref_grid',1.0,validScalarPosNum)
            parse(p,varargin{:});
            this.rho_ref     = p.Results.rho_ref;
            this.a_ref       = p.Results.a_ref;
            this.gamma       = p.Results.gamma;
            this.p_inf       = p.Results.p_inf / (this.rho_ref * this.a_ref^2);
            this.v_inf       = p.Results.v_inf / this.a_ref;
            this.l_ref       = p.Results.l_ref;
            this.l_ref_grid  = p.Results.l_ref_grid;
            this.alpha       = deg2rad( p.Results.alpha );
            this.circulation = p.Results.circulation;

            L = this.l_ref / this.l_ref_grid;
            this.a           = p.Results.a / L;
            this.b           = p.Results.b / L;
            this.R           = ( this.a + this.b )/2;
            % this.l           = (this.a^2 -this.b^2)/4;
            this.l = sqrt( (this.a^2 - this.b^2) / 4 );
        end
        function z = z_from_zeta(this,zeta)
            z = potential_ellipse.joukowsky_transform(this.l,zeta);
        end
        function dzdzeta = diff_z_from_zeta(this,zeta)
            dzdzeta = potential_ellipse.joukowsky_transform_derivative(this.l,zeta);
        end
        function zeta = zeta_from_z(this,z)
            zeta = potential_ellipse. ...
                   inverse_joukowsky_transform( this.l, z );
        end
        function F = complex_potential_from_zeta(this,zeta)
            F = potential_ellipse.cylinder_potential( this.v_inf,       ...
                                                      this.R,           ...
                                                      this.alpha,       ...
                                                      this.circulation, ...
                                                      zeta );
        end
        function w = complex_velocity_from_zeta(this,zeta)
            dzdzeta = this.diff_z_from_zeta(zeta);
            w = potential_ellipse.cylinder_velocity( this.v_inf,       ...
                                                     this.R,           ...
                                                     this.alpha,       ...
                                                     this.circulation, ...
                                                     zeta )./dzdzeta;
        end
        function F = complex_potential( this, z )
            F = this.complex_potential_from_zeta( this.zeta_from_z(z) );
        end
        function w = complex_velocity( this, z )
            w = this.complex_velocity_from_zeta( this.zeta_from_z(z) );
        end
        function rho = density( this, x, ~ )
            rho = this.rho_inf*ones(size(x));
        end
        function u = x_velocity( this, x, y )
            u = real( this.complex_velocity(x+1i*y) );
        end
        function v = y_velocity( this, x, y )
            v = -imag( this.complex_velocity(x+1i*y) );
        end
        function p = pressure( this, x, y )
            w2 = abs(this.complex_velocity(x+1i*y)).^2;
            p  = this.p_inf + 0.5*this.rho_inf*( this.v_inf^2 - w2 );
        end
        function h = plot_ellipse(this)
            h = fplot(@(theta)this.a*cos(theta),@(theta)this.b*sin(theta),[0,2*pi]);
        end
        function GRID = generate_elliptic_grid(this,xi_mult,n_eta,n_xi)
            c = sqrt( this.a^2 + this.b^2 );
            xi_0 = 0.5*log((this.a+this.b)/(this.a-this.b));
            % xi_0 = 0.5*log(abs((this.a+this.b)/(this.a-this.b)));
            xi_max = xi_mult*xi_0;
            xi  = linspace(xi_0,xi_max,n_xi);
            eta = linspace(0,2*pi,n_eta);
            [ETA,XI] = ndgrid(eta,xi);
            ZETA = XI + 1i*ETA;
            Z    = c*cosh(ZETA);
            GRID = struct();
            GRID.x = real(Z);
            GRID.y = imag(Z);
            GRID.imax = n_eta;
            GRID.jmax = n_xi;
        end
        function GRID = generate_mapped_grid(this,boundary_distance,imax,jmax)
            

            beta  = jmax - 1;
            eta_a = (this.a+boundary_distance)/this.a;
            alpha_a = eta_a^(1/beta);

            eta_b = (this.b+boundary_distance)/this.b;
            alpha_b = eta_b^(1/beta);

            t = linspace(0,1,jmax);
            theta = linspace(0,2*pi,imax).';

            a_ = this.a*alpha_a.^(beta*t);
            b_ = this.b*alpha_b.^(beta*t);

            GRID = struct();
            GRID.x = a_.*cos(theta);
            GRID.y = b_.*sin(theta);
            GRID.imax = imax;
            GRID.jmax = jmax;
        end
        function GRID = extruded_grid(this,imax,jmax)
            GRID = struct();
            GRID.imax = imax;
            GRID.jmax = jmax;
            GRID.x = zeros(imax,jmax);
            GRID.y = zeros(imax,jmax);

            theta = linspace(2*pi,0,imax);
            GRID.x(:,1) = this.a*cos(theta);
            GRID.y(:,1) = this.b*sin(theta);
            i21  = (imax-1)/2+1;
            i41  = (imax-1)/4+1;
            i34 = 3*(imax-1)/4+1;
            for j = 2:jmax
                [GRID.x(:,j),GRID.y(:,j)] = potential_ellipse.extrude_surface_pts(GRID.x(:,j-1),GRID.y(:,j-1),true);
                GRID.x(1,j)   = (GRID.x(1,j)+GRID.x(end,j))/2;
                GRID.x(end,j) =  GRID.x(1,j);
                GRID.x(i2,j)  = -GRID.x(1,j);

                GRID.y(1,j)    = 0;
                GRID.y(i2,j)   = 0;
                GRID.y(end,j)  = 0;
            end
        end
    end
    methods (Static)
        function F = cylinder_potential(v_inf,R,alpha,circulation,z)
            F = v_inf * ( z*exp(-1i*alpha) + R^2*exp(1i*alpha)./z ) ...
                      + (0.5*1i*circulation/pi)*log(z);
        end
        function w = cylinder_velocity(v_inf,R,alpha,circulation,z)
            w = v_inf * ( exp(-1i*alpha) - R^2*exp(1i*alpha)./z.^2 ) ...
                + (0.5*1i*circulation/pi)./z;
        end
        function z = joukowsky_transform(l,zeta)
            tol  = 1.0e-12;
            z    = zeros(size(zeta));
            mask = abs(zeta)>tol;
            z(mask)  = zeta(mask) + l^2./zeta(mask);
        end
        function dzdzeta = joukowsky_transform_derivative(l,zeta)
            tol     = 1.0e-12;
            dzdzeta = ones(size(zeta));
            mask    = abs(zeta)>tol;
            dzdzeta(mask) = 1 - l^2./zeta(mask).^2;
        end
        function zeta = inverse_joukowsky_transform(l,z)
            tol  = 1.0e-12;
            if l>tol
                zp   = sqrt(z+2*l);
                zm   = sqrt(z-2*l);
                zeta = l*( (zp+zm)./(zp-zm) );
            else
                zeta = z;
            end
        end
        function [x,y,h] = extrude_surface_pts(xs,ys,periodic,h)
            % if (periodic)
            %     s = sign(polygon_area(xs(1:end-1),ys(1:end-1)));
            % else
                % s = sign(polygon_area(xs,ys));
            % end
            % if s == 0
            %     error('could not determine normal vector')
            % end
            s = -1;
            if (periodic)
                dx = gradient([xs(end-1);xs(:);xs(2)]);
                dx = dx(2:end-1);
                dy = gradient([ys(end-1);ys(:);ys(2)]);
                dy = dy(2:end-1);
            else
                dx = gradient(xs(:));
                dy = gradient(ys(:));
            end
            mag = sqrt( dx.^2 + dy.^2 );

            % fix for singular points
            mag(mag<eps(1)) = 1;

            % outward facing normal
            n1 =  s*dy./mag;
            n2 = -s*dx./mag;
            if nargin<4
                % h = min(mag);
                h = mag;
                % h1 = min(mag);
                % h = h1*log((h/h1)*(exp(1)-1) + 1);
                % h = h1*sqrt(h/h1);
            end
            % add offset
            x =xs + n1.*h(:);
            y =ys + n2.*h(:);
            function A = polygon_area(x,y)
                x1 = [x(:);x(1)];
                y1 = [y(:);y(1)];
                A = 0.5*sum( x1(1:end-1).*y1(2:end) - x1(2:end).*y1(1:end-1) );
            end
        end
    end
end