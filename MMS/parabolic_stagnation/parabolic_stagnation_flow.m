classdef parabolic_stagnation_flow
    properties
        r_LE(1,1)        double = 1.0
        eta_0(1,1)       double = 1.0
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
        function this = parabolic_stagnation_flow(varargin)
            validScalarNum        = @(x) isnumeric(x) && isscalar(x);
            validScalarNonNegNum  = @(x) validScalarNum(x) && (x >= 0);
            validScalarPosNum     = @(x) validScalarNum(x) && (x > 0);
            p = inputParser;
            p.addOptional('rLE',1.0,validScalarNonNegNum);
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
            L = this.l_ref / this.l_ref_grid;
            this.r_LE         = p.Results.rLE / L;
            this.eta_0        = sqrt( this.r_LE );
        end
        function z = z_from_zeta( this, zeta )
            z = 0.5*( zeta.^2 + this.r_LE );
        end
        function dz_dzeta = diff_z_from_zeta( ~, zeta )
            dz_dzeta = zeta;
        end
        function zeta = zeta_from_z( this, z )
            s = sign(imag(z));
            z2 = real(z) + 1i*abs(imag(z));
            zeta1 = sqrt( 2*z2 - this.r_LE );
            zeta = s.*real(zeta1) + 1i*imag(zeta1);
        end
        function dzeta_dz = diff_zeta_from_z( this, z )
            dzeta_dz = 1./sqrt( 2*z - this.r_LE );
        end
        function r = r_from_theta( this, theta )
            r = this.r_LE./(1-cos(theta));
        end
        function x = x_from_theta(this,theta)
            x = this.r_from_theta(theta).*cos(theta) + 0.5*this.r_LE;
        end
        function y = y_from_theta(this,theta)
            y = this.r_from_theta(theta).*sin(theta);
        end
        function theta = theta_from_x(this,x)
            a = (x/this.r_LE) - 0.5;
            theta = acos( a./(1+a) );
        end
        function y = y_from_x(this,x)
            y = sqrt( 2*this.r_LE*x );
        end
        function r = r_from_x(this,x)
            y = this.y_from_x(x);
            r = sqrt( x.^2 +  y.^2 );
        end
        function [x,y] = surface_coords_theta(this,theta)
            r = this.r_from_theta(theta);
            x =  r.*cos(theta) + 0.5*this.r_LE;
            y =  r.*sin(theta);
        end
%% Plot
        function [h1] = plot_surface_theta(this,theta_min,varargin)
            h1 = fplot( @(theta)this.x_from_theta(theta), ...
                        @(theta)this.y_from_theta(theta), ...
                        [theta_min,2*pi-theta_min], varargin{:} );
            set(gca,'DataAspectRatio',[1 1 1])
        end
        function [h1] = plot_surface_x(this,x_max,varargin)
            h1 = fplot( @(t)abs(t),@(t)sign(t).*this.y_from_x(abs(t)), ...
                [-x_max,x_max], varargin{:} );
            set(gca,'DataAspectRatio',[1 1 1])
        end
        function w = complex_velocity(this,zeta)
            w = this.v_inf * ( 1 - 1i*this.eta_0./zeta );
        end
        function u = x_velocity(this,x,y)
            u =  real( this.complex_velocity(this.zeta_from_z(x+1i*y)));
        end
        function v = y_velocity(this,x,y)
            v = -imag( this.complex_velocity(this.zeta_from_z(x+1i*y)));
        end
        function p = pressure(this,x,y)
            w  = this.complex_velocity(this.zeta_from_z(x+1i*y));
            p  = this.p_inf + (1/2)*this.rho_inf*( this.v_inf^2 - abs(w).^2 );
        end
        function ut = surface_tangent_velocity(this,x)
            ut = this.v_inf * sqrt( 2*x ./ ( 2*x + this.r_LE) );
        end
        % function ut = surface_tangent_velocity_2(this,x)
        %     y = this.y_from_x(x);
        %     w = this.complex_velocity( this.zeta_from_z( x + 1i*y ) );
        %     ut = abs(w);
        % end
        function p = surface_pressure_x(this,x)
            ut = this.surface_tangent_velocity(x);
            p  = this.p_inf + (1/2)*this.rho_inf*( this.v_inf^2 - ut.^2 );
        end
        function p = surface_pressure_theta(this,theta)
            x  = this.x_from_theta(theta);
            ut = this.surface_tangent_velocity(x);
            p  = this.p_inf + (1/2)*this.rho_inf*( this.v_inf^2 - ut.^2 );
        end
        function rho = density(this,x,~)
            rho = this.rho_inf*ones(size(x));
        end
        function s = arc_length(this,x)
            xi = sqrt(2*x);
            s = 0.5*(xi.*sqrt(xi.^2+this.eta_0^2) + this.eta_0^2*asinh(xi/this.eta_0));
        end
        function s = arc_length_segment(this,x0,x1)
            s = this.arc_length(x1) - this.arc_length(x0);
        end
        function x = arc_length_param_x(this,x_max,N,F)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            x  = zeros(N,1);
            t  = linspace(0,1,N).';
            x(1) = 0;
            x(N) = x_max; 
            L_total = this.arc_length(x_max);
            dL = (F(1)-F(0))/ L_total;
            for i = 2:N-1
                ftmp = @(xsi) F(t(i))-F(t(i-1)) - dL * this.arc_length_segment(x(i-1),xsi);
                x(i) = fzero( @(xsi)ftmp(xsi),[0,x_max],options);
            end
        end
        function GRID = extruded_grid(this,n_theta,n_r,stag_spacing,boundary_distance,AR)
            GRID = struct();
            GRID.imax = n_theta;
            GRID.jmax = n_r;
            GRID.x = zeros(n_theta,n_r);
            GRID.y = zeros(n_theta,n_r);
            L = this.arc_length(boundary_distance);
            h = stag_spacing/AR;
            d0  = stag_spacing/L;
            off = 0.05;
            f = hermite_blend_2_vinokur_one_sided(n_theta,0.5*d0,off,true);
            f2 = @(t) f(0.5 + t/2);
            N2 = (n_theta-1)/2+1;
            % generate surface points
            x = this.arc_length_param_x(boundary_distance,N2,f2);
            y = this.y_from_x(x);
            GRID.x(N2:n_theta,1) = x;
            GRID.x(N2-1:-1:1,1) = x(2:N2);
            GRID.y(N2:n_theta,1) = y;
            GRID.y(N2-1:-1:1,1) = -y(2:N2);           
            [~,alpha,~] = parabolic_stagnation_flow.geomspace( n_r, 0, boundary_distance, h );
            for j = 2:n_r
                h = alpha*h;
                [GRID.x(:,j),GRID.y(:,j)] = parabolic_stagnation_flow.extrude_surface_pts(GRID.x(:,j-1),GRID.y(:,j-1),h);
            end
            % check
            % max( abs( GRID.y(end:-1:N2,:)+GRID.y(1:N2,:) ),[], 'all')==0
            % max( abs( GRID.x(end:-1:N2,:)-GRID.x(1:N2,:) ),[], 'all')==0
        end
        function GRID = parabolic_grid(this,n_theta,n_r,stag_spacing,boundary_distance,AR)
            GRID = struct();
            GRID.imax = n_theta;
            GRID.jmax = n_r;
            GRID.x = zeros(n_theta,n_r);
            GRID.y = zeros(n_theta,n_r);
            L = this.arc_length(boundary_distance);
            h = stag_spacing/AR;
            d0  = stag_spacing/L;
            off = 0.05;
            f = hermite_blend_2_vinokur_one_sided(n_theta,0.5*d0,off,true);
            f2 = @(t) f(0.5 + t/2);
            N2 = (n_theta-1)/2+1;
            % generate surface points
            x = this.arc_length_param_x(boundary_distance,N2,f2);
            y = this.y_from_x(x);
            z = x + 1i*y;

            xi = real( this.zeta_from_z( z ) );
            eta0 = abs( this.zeta_from_z( z(1) ) );
            deta_0 = abs( this.zeta_from_z( z(1) - h) - this.zeta_from_z( z(1) ) );
            eta_1 = abs( this.zeta_from_z( z(1) - boundary_distance) - this.zeta_from_z( z(1) ) );
            [eta,~,~] = parabolic_stagnation_flow.geomspace( n_r, eta0, eta_1, deta_0 );

            [XI,ETA] = ndgrid(xi,eta);
            ZETA = XI + 1i * ETA;

            % zeta1 = xi(1)   + 1i*eta(1);
            % zeta2 = xi(end) + 1i*eta(1);
            % zeta3 = this.zeta_from_z( z(end) + 1i*500);
            % zeta4 = xi(1)   + 1i*eta(end);
            % ZETA = parabolic_stagnation_flow.bilinear_remap_rectangle_cmplx(zeta1,zeta2,zeta3,zeta4,xi,eta);
            
            Z = this.z_from_zeta(ZETA);
            X = real(Z);
            Y = imag(Z);
            
            GRID.x(N2:n_theta,:) = X;
            GRID.x(N2-1:-1:1,:) = X(2:N2,:);
            GRID.y(N2:n_theta,:) = Y;
            GRID.y(N2-1:-1:1,:) = -Y(2:N2,:);
        end
    end
    methods (Static)
        function ZETA = bilinear_remap_rectangle_cmplx(zeta1,zeta2,zeta3,zeta4,xi0,eta0)
            xi  = 2*(xi0-xi0(1))/(xi0(end)-xi0(1)) - 1;
            eta = 2*(eta0-eta0(1))/(eta0(end)-eta0(1)) - 1;
            [XI,ETA] = ndgrid(xi,eta);
            % shape functions
            P1 = @(xi,eta) 0.25*(1-xi).*(1-eta);
            P2 = @(xi,eta) 0.25*(1+xi).*(1-eta);
            P3 = @(xi,eta) 0.25*(1+xi).*(1+eta);
            P4 = @(xi,eta) 0.25*(1-xi).*(1+eta);

            % corner points
            ZETA = zeta1*P1(XI,ETA) ...
                 + zeta2*P2(XI,ETA) ...
                 + zeta3*P3(XI,ETA) ...
                 + zeta4*P4(XI,ETA);
        end
        function [x,y,h] = extrude_surface_pts(xs,ys,h)
            s = sign(polygon_area(xs,ys));
            if s == 0
                error('could not determine normal vector')
            end
            dx = gradient(xs(:));
            dy = gradient(ys(:));
            mag = sqrt( dx.^2 + dy.^2 );
            
            % outward facing normal
            n1 =  s*dy./mag;
            n2 = -s*dx./mag;
            if nargin<3
                h = min(mag);
                % h = mag;
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
        function [t,tc,L] = reparam_curve_xy(f1,f2,N,tmin,tmax,mu)
            % creates a set of points parameterized in (approximate) arc length space
            % of f1(tmin,tmax) with spacing from f2(0,1)
            t0 = linspace(0,1,N);
            [x,y] = f1( tmin + (tmax-tmin)*t0 );
            points = [x(:).';y(:).'];
            tc = [ 0; cumsum( sqrt( sum( abs(points(:,2:end) - points(:,1:end-1)).^2, 1).^mu ) ).'];
            L = tc(end);
            tc = tc/L;
            t = interp1(tc,t0,f2(t0),"spline");
            t = tmin + (tmax-tmin)*t;
        end
        function [x,r,dx1] = geomspace( N, xmin, xmax, dx0 )
            % geometric spacing -> each subsequent interval is r times the length of
            % previous
            % e.g. r = 1.1 gives a 10% increase in delta x for adjacent nodes
            
            S = xmax - xmin;
            
            if S/(N-2) > dx0
                r0 = 1.01;
            else
                r0 = 0.99;
            end
            
            fun = @(r) S/dx0 - sum(r.^(0:N-2));
            options = optimset('FunValCheck','on');
            r = fzero(fun,r0,options);
            
            x = zeros(1,N);
            x(1) = xmin;
            for i = 2:N-1
                x(i) = x(i-1) + (r^(i-2)*dx0);
            end
            x(N) = xmax;
            dx1 = x(N)-x(N-1);
        end
    end
end