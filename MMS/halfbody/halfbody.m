classdef halfbody
    properties
        l(1,1)         double  = 1.0
        sigx2pi(1,1)   double  = 0.0
        xstag(1,1)     double  = 0.0
        ymax(1,1)      double  = 0.0
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
            % theta = wrapTo2Pi(theta);
            y = (this.sigx2pi/this.vinf)*(pi-theta) * this.l;
        end
        function r = r_from_theta(this,theta)
            % theta = wrapTo2Pi(theta);
            tol = 1e-5;
            r0 = -this.xstag;
            r1 = this.y_from_theta(pi-tol)./sin(pi-tol);

            theta(theta>pi) = 2*pi-theta(theta>pi);
            r = r0*ones(size(theta));
            mask = abs(pi-theta)>tol;
            
            r(mask) = this.y_from_theta(theta(mask))./sin(theta(mask));
            r(~mask) = r0 + (r1-r0)*(pi-theta(~mask))/tol;
        end
        function x = x_from_theta(this,theta)
            x = this.r_from_theta(theta).*cos(theta);
        end
        function [theta] = theta_from_x(this,x)
            theta_0 = asin(this.ymax/x);
            theta = fzero(@(theta)this.x_from_theta(theta)-x,theta_0);
        end
        function [x,y] = surface_coords(this,theta)
            r = this.r_from_theta(theta);
            x = r.*cos(theta);
            y = r.*sin(theta);
        end
%% Plot
        function [h1] = plot_halfbody(this,theta_min,varargin)
            h1 = fplot( @(theta)this.x_from_theta(theta), ...
                        @(theta)this.y_from_theta(theta), ...
                        [theta_min,2*pi-theta_min], varargin{:} );
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
        function rho = density(this,x,~)
            rho = this.rhoinf*ones(size(x));
        end
        function GRID = extruded_grid(this,n_theta,n_r,boundary_distance,mu,AR)
            % estimate theta_min
            % theta_min1 = asin(this.ymax/boundary_distance);
            theta_min = this.theta_from_x(boundary_distance);

            theta_max = 2*pi-theta_min;
            GRID = struct();
            GRID.imax = n_theta;
            GRID.jmax = n_r;
            GRID.x = zeros(n_theta,n_r);
            GRID.y = zeros(n_theta,n_r);
            % theta = linspace(theta_min,theta_max,n_theta);
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