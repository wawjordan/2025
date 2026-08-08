classdef halfbody
    properties
        l(1,1)         double  = 1.0
        sigx2pi(1,1)   double  = 0.0
        xstag(1,1)     double  = 0.0
        ymax(1,1)      double  = 0.0
        theta_cp0(1,1) double = 1.97603146838258 % 0 = tan(t) - 2*(t-pi)
        theta_cpmin(1,1) double = 1.09880571085138 % 0 = 1 + (2*(pi-t)^2-1)*cos(2*t) + 2*(pi-t)*sin(2*t)
        rhoinf(1,1)    double  = 1.0
        vinf(1,1)      double  = 1.0
        pinf(1,1)      double  = 1.0
        gamma(1,1)     double  = 1.4
    end  
    methods
%% Constructor
        function this = halfbody(xstag,vinf,l)
            % we'll go with the convention of velocity always moving to the right
            this.xstag = -abs(xstag);
            if ( nargin > 1)
                this.vinf  = vinf;
            end
            if ( nargin > 2)
                this.l = l;
            end
            this.sigx2pi = -this.vinf*this.xstag;
            this.ymax  = this.sigx2pi*pi/this.vinf;
        end
        function y = y_from_theta(this,theta)
            y = (this.sigx2pi/this.vinf)*(pi-theta) * this.l;
        end
        function r = r_from_theta(this,theta)
            r = this.l*(this.sigx2pi/this.vinf)*(pi-theta)./sin(theta);
        end
        function r = r_from_theta_safe(this,theta)
            tol = 1e-5;
            r0 = -this.xstag;
            r1 = this.r_from_theta(pi-tol);
            theta(theta>pi) = 2*pi-theta(theta>pi);
            r = r0*ones(size(theta));
            mask = abs(pi-theta)>tol;
            r(mask) = this.r_from_theta(theta(mask));
            r(~mask) = r0 + (r1-r0)*(pi-theta(~mask))/tol;
        end
        function drdtheta = diff_r_from_theta(this,theta)
            c = this.l*(this.sigx2pi/this.vinf);
            drdtheta = -c * ( (pi-theta).*cos(theta) + sin(theta) )./(sin(theta).^2);
        end
        function d2rdtheta2 = diff2_r_from_theta(this,theta)
            c = this.l*(this.sigx2pi/this.vinf);
            ct = cos(theta);
            st = sin(theta);
            d2rdtheta2 = 2*c*ct./st.^2 - c*(theta-pi)./st - 2*c*ct.^2.*(theta-pi)./st.^3;
        end
        function d2rdtheta2 = diff2_r_from_theta_safe(this,theta)
            tol = 1e-5;
            theta(theta>pi) = 2*pi-theta(theta>pi);
            d2rdtheta2 = zeros(size(theta));
            mask = abs(pi-theta)>tol;
            d2rdtheta2(mask) = this.diff2_r_from_theta(theta(mask));
            d2rdtheta2(~mask) = 0;
        end
        function rho = curvature(this,theta)
            r   = this.r_from_theta_safe(theta);
            dr  = this.diff_r_from_theta_safe(theta);
            ddr = this.diff2_r_from_theta_safe(theta);
            rho = abs( r.^2 + 2*dr.^2 - r.*ddr )./(r.^2 + dr.^2).^(3/2);
        end
        function drdtheta = diff_r_from_theta_safe(this,theta)
            tol = 1e-5;
            r0 = -this.xstag;
            r1 = this.r_from_theta(pi-tol);
            theta(theta>pi) = 2*pi-theta(theta>pi);
            drdtheta = zeros(size(theta));
            mask = abs(pi-theta)>tol;
            drdtheta(mask) = this.diff_r_from_theta(theta(mask));
            drdtheta(~mask) = -(r1-r0)/tol;
        end
        function ds = diff_arc_length(this,theta)
            r = this.r_from_theta_safe(theta);
            drdtheta = this.diff_r_from_theta_safe(theta);
            ds = sqrt( r.^2 + drdtheta.^2 );
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function S = arc_length(this,t0,t1)
            S = integral(@(t)this.diff_arc_length(t),t0,t1,"AbsTol",1e-15,'RelTol',1e-14);
        end
        function ts = arc_length_param(this,t,F)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            N = numel(t);
            L_total = this.arc_length(t(1),t(N));
            dL = ( t(N) - t(1) ) / L_total;
            ts = zeros(N,1);
            ts(1) = t(1);
            ts(N) = t(N);
            for i = 2:N-1
                ftmp = @(tsi) F(t(i))-F(t(i-1)) - dL * this.arc_length(ts(i-1),tsi);
                ts(i) = fzero( @(tsi)ftmp(tsi),ts(i-1),options);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = x_from_theta(this,theta)
            x = this.r_from_theta_safe(theta).*cos(theta);
        end
        function [theta] = theta_from_x(this,x)
            theta_0 = asin(this.ymax/x);
            theta = fzero(@(theta)this.x_from_theta(theta)-x,theta_0);
        end
        function [x,y] = surface_coords(this,theta)
            r = this.r_from_theta_safe(theta);
            x = r.*cos(theta);
            y = r.*sin(theta);
        end
        function [n1,n2] = unit_normal(this,theta)
            nr =  this.r_from_theta_safe(theta);
            nt = -this.diff_r_from_theta_safe(theta);
            nt(theta>pi) = -nt(theta>pi);
            n1 = cos(theta).*nr - sin(theta).*nt;
            n2 = sin(theta).*nr + cos(theta).*nt;
            mag = sqrt(n1.^2+n2.^2);
            n1 = n1./mag;
            n2 = n2./mag;
        end
        function [Fx,Fy] = forces(this,t0,t1)
            Fx = this.surface_integral(@(theta)px_force(this,theta),t0,t1);
            Fy = this.surface_integral(@(theta)py_force(this,theta),t0,t1);
            function px = px_force(this,theta)
                [px,~] = this.pressure_force(theta);
            end
            function py = py_force(this,theta)
                [~,py] = this.pressure_force(theta);
            end
        end
        function val = surface_integral(this,fun,t0,t1)
            val = integral(@(theta)fun(theta) ...
                           .*this.diff_arc_length(theta),...
                           t0,t1,"AbsTol",1e-14,'RelTol',1e-12);
        end
        function [pnx,pny] = pressure_force(this,theta)
            [nx,ny] = this.unit_normal(theta);
            p = (this.pinf-this.surface_pressure(theta));
            pnx = p.*nx;
            pny = p.*ny;
        end
%% Plot
        function [h1] = plot_halfbody(this,theta_min,varargin)
            h1 = fplot( @(theta)this.x_from_theta(theta), ...
                        @(theta)this.y_from_theta(theta), ...
                        [theta_min,2*pi-theta_min], varargin{:} );
            set(gca,'DataAspectRatio',[1 1 1])
        end
        function [h1] = plot_halfbody_normals(this,theta,varargin)
            [n1,n2] = this.unit_normal(theta);
            [x,y]   = this.surface_coords(theta);
            h1 = quiver(x,y,n1,n2,varargin{:});
            set(gca,'DataAspectRatio',[1 1 1])
        end
        function u = x_velocity(this,x,y)
            u = this.vinf + this.sigx2pi*x./(x.^2+y.^2);
        end
        function v = y_velocity(this,x,y)
            v = this.sigx2pi*y./(x.^2+y.^2);
        end
        function p = pressure(this,x,y)
            u = this.x_velocity(x,y);
            v = this.y_velocity(x,y);
            p  = this.pinf + (1/2)*this.rhoinf*( this.vinf^2 - (u.^2 + v.^2) );
        end
        function p = surface_pressure(this,theta)
            r = this.r_from_theta_safe(theta);
            u =  this.vinf + this.sigx2pi*cos(theta)./r;
            v = -this.sigx2pi*sin(theta)./r;
            p  = this.pinf + (1/2)*this.rhoinf*( this.vinf^2 - (u.^2 + v.^2) );
        end
        function rho = density(this,x,~)
            rho = this.rhoinf*ones(size(x));
        end
        function GRID = extruded_grid1(this,n_theta,n_r,boundary_distance,mu,AR)
            theta_min = this.theta_from_x(boundary_distance);
            theta_max = 2*pi-theta_min;
            GRID = struct();
            GRID.imax = n_theta;
            GRID.jmax = n_r;
            GRID.x = zeros(n_theta,n_r);
            GRID.y = zeros(n_theta,n_r);
            f = @(t) t;
            [theta,~,~] = halfbody.reparam_curve_xy(@(theta)this.surface_coords(theta),f,n_theta,theta_min,theta_max,mu);
            % generate surface points
            [GRID.x(:,1),GRID.y(:,1)] = this.surface_coords(theta);
            [~,~,h] = halfbody.extrude_surface_pts(GRID.x(:,1),GRID.y(:,1));
            h = h/AR;
            [~,alpha,~] = halfbody.geomspace( n_r, abs(this.xstag), boundary_distance, h );
            for j = 2:n_r
                h = alpha*h;
                [GRID.x(:,j),GRID.y(:,j)] = halfbody.extrude_surface_pts(GRID.x(:,j-1),GRID.y(:,j-1),h);
            end
        end
        function GRID = extruded_grid2(this,n_theta,n_r,stag_spacing,boundary_distance,AR)
            theta_min = this.theta_from_x(boundary_distance);
            theta_max = 2*pi-theta_min;
            GRID = struct();
            GRID.imax = n_theta;
            GRID.jmax = n_r;
            GRID.x = zeros(n_theta,n_r);
            GRID.y = zeros(n_theta,n_r);
            % f = @(t) t;
            L = this.arc_length(theta_min,theta_max);
            h = stag_spacing/AR;
            d0  = stag_spacing/L;
            off = 0.05;%10*d0;
            f  = hermite_blend_2_vinokur_one_sided(n_theta,d0,off,true);
            f = @(theta) theta_min + (theta_max-theta_min)*f((theta-theta_min)/(theta_max-theta_min));
            N2 = (n_theta-1)/2+1;
            theta0 = linspace(theta_min,pi,N2);
            theta = zeros(n_theta,1);
            theta(1:N2) = this.arc_length_param(theta0,f);
            theta(N2+1:end) = 2*pi-theta(N2-1:-1:1);
            % generate surface points
            [GRID.x(:,1),GRID.y(:,1)] = this.surface_coords(theta);
            [~,alpha,~] = halfbody.geomspace( n_r, abs(this.xstag), boundary_distance, h );
            for j = 2:n_r
                h = alpha*h;
                [GRID.x(:,j),GRID.y(:,j)] = halfbody.extrude_surface_pts(GRID.x(:,j-1),GRID.y(:,j-1),h);
            end
        end
        % function [mu,t,tc,L] = get_optimal_mu(f1,N,theta_min)
        %     [~,tc,L] = halfbody.reparam_curve_xy(@(theta)this.surface_coords(theta),f,n_theta,theta_min,theta_max,mu);
        %     min(diff(tc))*L
        % end
    end
    methods (Static)
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