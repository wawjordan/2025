classdef kt_airfoil
    properties
        l(1,1)         double  = 1.0
        epsilon(1,1)   double  = 0.0
        kappa(1,1)     double  = 0.0
        tau(1,1)       double  = 0.0
        rhoinf(1,1)    double  = 1.0
        vinf(1,1)      double  = 1.0
        pinf(1,1)      double  = 1.0
        gamma(1,1)     double  = 1.4
        chord(1,1)     double  = 1.0
        thetaLE(1,1)   double  = 0.0
        thetaCmax(1,1) double  = 0.0
        thetaSP(1,1)   double  = 0.0;
        thetaTE(1,1)   double  = 0.0
        zLE(1,1)       double  = 0.0 + 1i*0.0
        zTE(1,1)       double  = 0.0 + 1i*0.0
        xLE(1,1)       double  = 0.0
        yLE(1,1)       double  = 0.0
        xTE(1,1)       double  = 0.0
        yTE(1,1)       double  = 0.0
        n(1,1)         double  = 2.0
        alpha(1,1)     double  = 0.0
        beta(1,1)      double  = 0.0
        a(1,1)         double  = 1.0
        mu(1,1)        double  = 0.0 + 1i*0.0
        ref_z(:,1)     double  = 0.0 + 1i*0.0
        ref_t(:,1)     double  = 0.0 + 1i*0.0
        n_ref(1,1)     double  = 4096
        bounds(4,1)    double  = 0.0
    end  
    methods
%% Constructor
        function this = kt_airfoil(epsilon,kappa,tau,n_ref)
            this.epsilon = epsilon;
            this.kappa   = kappa;
            this.tau     = tau;
            this.n       = 2 - this.tau/pi;
            this.a       = this.l*sqrt( (1+this.epsilon)^2 + this.kappa^2 );
            this.mu      = this.l*(-this.epsilon + this.kappa*1i);
            this.beta    = asin(this.l*this.kappa/this.a);
            this         = this.set_airfoil_geometry();
            if ( nargin > 3)
                this.n_ref = n_ref;
            end
        end
%% Set airfoil geometry
        function this = set_airfoil_geometry(this)
            % angle at leading edge
            this.thetaLE = fminbnd(@(theta)real(this.z_from_theta(theta)),0,2*pi);
            this.thetaTE = 0;
            % NOTE: These are not scaled for unit chord
            this.zLE = this.z_from_theta(this.thetaLE);
            this.zTE = this.z_from_theta(this.thetaTE);
            this.chord = ( abs(this.zLE) + abs(this.zTE) );
            this.xLE   = real(this.zLE);
            this.yLE   = imag(this.zLE);
            this.xTE   = real(this.zTE);
            this.yTE   = imag(this.zTE);

            % reference theta and z for use as a look-up table
            this.ref_t = linspace(0,2*pi,this.n_ref);
            this.ref_z = this.z_from_theta(this.ref_t);

            % bounding box for airfoil
            this.bounds(1) = min( real(this.ref_z) ) - 0.1;
            this.bounds(2) = max( real(this.ref_z) ) + 0.1;
            this.bounds(3) = min( imag(this.ref_z) ) - 0.1;
            this.bounds(4) = max( imag(this.ref_z) ) + 0.1;

            % maximum curvature (tolerance to keep away from TE singularity)
            tol = 0.01;
            c_fun = @(theta)-this.airfoil_curvature2(theta);
            [this.thetaCmax,~] = fminbnd(c_fun,0+tol,2*pi-tol);
            % leading edge stagnation point (changed in set_alpha)
            this.thetaSP = pi;
        end
%% Set angle of attack
        function this = set_alpha(this,alpha)
            % angle of attack (specified in degrees)
            this.alpha = deg2rad(alpha);
            % calculate leading edge stagnation point
            tol = 0.01;
            this.thetaSP = fminbnd( @(theta) (abs(this.w_cylinder(this.zeta_from_theta(theta)))).^2, 0+tol, 2*pi-tol );
        end
%% Analytic Force Coefficients
        function val = CL(this)
            % analytic lift coefficient
            val = 8*pi*(this.a/this.chord)*sin(this.alpha+this.beta);
        end
        function val = CD(~)
            % analytic drag coefficient
            val = 0;
        end
        function val = CX(this)
            % analytic axial force coefficient
            val = -this.CL * sin(this.alpha);
        end
        function val = CY(this)
            % analytic normal force coefficient
            val = this.CL * cos(this.alpha);
        end
%% Plot airfoil
        function [h1,h2] = plot_airfoil(this,scale)
            if (scale)
                zfun = @(theta) ( this.z_from_theta(theta) - this.zLE ) ...
                                ./this.chord;
            else
                zfun = @(theta) this.z_from_theta(theta);
            end
            h1 = fplot( @(theta) real(zfun(theta)), ...
                        @(theta) imag(zfun(theta)), [0,this.thetaLE] );
            hold on;
            h2 = fplot( @(theta) real(zfun(theta)), ...
                        @(theta) imag(zfun(theta)), [this.thetaLE,2*pi] );
            set(gca,'DataAspectRatio',[1 1 1])
        end
        function [h] = plot_osculating_circles(this,theta,n,scale)
            t = linspace(0,2*pi,n);
            h = struct();
            [r,x,y] = this.get_osculating_circle(theta,scale);
            for i = 1:numel(theta)
                h.p(i) = plot(x(i)+r(i)*cos(t),y(i)+r(i)*sin(t));
            end
        end
        function [r,x,y] = get_osculating_circle(this,theta,scale)
            K = this.airfoil_curvature2(theta);
            r = 1./K;
            [N1,N2,~] = this.unit_normal(theta);
            if scale
                z = ( this.z_from_theta(theta) - this.zLE )./this.chord;
            else
                z = this.z_from_theta(theta);
            end
            x = real(z(:)) - r(:).*N1(:);
            y = imag(z(:)) - r(:).*N2(:);
        end
        function [h1,h2] = plot_cp(this)
            tol = 1e-12; % tolerance to avoid issues at trailing edge
            xfun = @(theta) real( this.z_from_theta(theta) - this.zLE ) ...
                            ./ this.chord;
            pfun = @(theta) ( this.airfoil_pressure( ...
                                   this.zeta_from_theta(theta) ) ...
                                   - this.pinf ) ...
                            ./ (0.5*this.rhoinf*this.vinf^2);
            h1 = fplot( @(theta) xfun(theta), ...
                        @(theta) pfun(theta), [0+tol,this.thetaLE] );
            hold on;
            h2 = fplot( @(theta) xfun(theta), ...
                        @(theta) pfun(theta), [this.thetaLE,2*pi-tol] );
        end
%% Scale and Shift
%% ( z --> zs )
        function z = scale_to_unit_chord(this,z)
            % scale and shift the coordinate z to be on an airfoil with unit
            % chord and LE at the origin
            z = (z - this.zLE)/this.chord;
        end
%% ( zs --> z )
        function z = scale_from_unit_chord(this,z)
            % scale and shift the coordinate z from an airfoil with unit
            % chord and LE at the origin to the original coordinates in the
            % z = z(zeta(theta)) plane
            z = z*this.chord + this.zLE;
        end
%% Mappings
%% zeta from theta ( theta --> zeta ) (on circle in complex plane)
        function zeta = zeta_from_theta(this,theta)
            % map to cylinder in complex plane
            zeta = this.a*exp(1i*(theta-this.beta)) + this.mu;
        end
%% Karman Trefftz Map  ( zeta --> z )
        function z = z_from_zeta(this,zeta)
            zeta_p = (zeta+this.l).^this.n;
            zeta_m = (zeta-this.l).^this.n;
            z = this.n*this.l*( zeta_p + zeta_m )./( zeta_p - zeta_m );
        end
%% Scaled and Shifted KT Map ( zeta --> z --> zs )
        function zs = zs_from_zeta(this,zeta)
            zs = this.scale_to_unit_chord( this.z_from_zeta(zeta) );
        end
%% Inverse Karman Trefftz Map  ( z --> zeta )
        function zeta = zeta_from_z(this,z)
            % NOTE: this isn't always 1 to 1
            z_p = (z+this.n*this.l);
            z_m = (z-this.n*this.l);
            xn = 1/this.n;
            zeta = -this.l*( (z_m./z_p).^xn + 1) ./ ( (z_m./z_p).^xn - 1);
        end
%% Scaled and Shifted Inverse KT Map  ( zs --> z --> zeta )
        function zeta = zeta_from_zs(this,zs)
            zeta = this.zeta_from_z( this.scale_from_unit_chord(zs) );
        end
%% Numeric Inverse Karman Trefftz Map ( z --> zeta )
        function zeta = zeta_from_z_num(this,z)
            % NOTE: this isn't always robust
            zeta = zeros(size(z));
            mask = in_bounds(this,z);
            zeta(~mask) = this.zeta_from_z(z(~mask));
            zeta(mask)  = arrayfun( @(z) zeta_lu_from_theta(this,z), z(mask) );
            function in = in_bounds(this,z)
                in = ( ( real(z)>this.bounds(1) )|( real(z)<this.bounds(2) ) )...
                    &( ( imag(z)>this.bounds(3) )|( imag(z)<this.bounds(4) ) );
            end
            function zeta = zeta_lu_from_theta(this,z)
               % query look-up table for initial guess
               [~,idx] = min(abs(this.ref_z - z));
               % solve using initial guess
               options = optimset('TolFun',1e-15,'TolX',1e-17,'Display','none');
               zeta = fsolve( @(zeta) z - this.z_from_zeta(zeta), ...
                             this.zeta_from_theta( this.ref_t(idx) ), options );
            end
        end
%% Numeric Scaled and Shifted Inverse KT Map ( zs --> z --> zeta )
        function zeta = zeta_from_zs_num(this,zs)
            zeta = this.zeta_from_z_num( this.scale_from_unit_chord(zs) );
        end
%% z from theta ( theta --> zeta --> z ) (on surface of airfoil)
        function z = z_from_theta(this,theta)
            z = this.z_from_zeta( this.zeta_from_theta(theta) );
        end
%% zs from theta ( theta --> zeta --> z --> zs ) (on surface of scaled and shifted airfoil)
        function zs = zs_from_theta(this,theta)
            zs = this.scale_to_unit_chord( this.z_from_theta(theta) );
        end
%% zeta from theta and r ( [theta,r] --> zeta ) (exterior to circle in complex plane)
        function zeta = zeta_from_theta_r(this,theta,r)
            zeta = r.*exp(1i*(theta-this.beta)) + this.mu;
        end
%% z from theta and r ( [theta,r] --> zeta --> z ) (exterior to airfoil)
        function z = z_from_theta_r(this,theta,r)
            z = this.z_from_zeta( this.zeta_from_theta_r(theta,r));
        end
%% z from theta and r ( [theta,r] --> zeta --> z --> zs ) (exterior to scaled and shifted airfoil)
        function z = zs_from_theta_r(this,theta,r)
            z = this.scale_to_unit_chord( this.z_from_theta_r(theta,r) );
        end
%% theta from zeta ( zeta --> [theta,r] )
        function [theta,r] = theta_from_zeta(this,zeta)
            zeta_ = zeta-this.mu;
            r = abs(zeta-this.mu);
            % theta = real(-1i*log(zeta_./r)+this.beta);
            theta = atan2(imag(zeta_),real(zeta_))+this.beta;
            theta = mod(theta,2*pi);
        end
%% theta from z ( z --> zeta --> [theta,r] )
        function [theta,r] = theta_from_z(this,z)
            [theta,r] = this.theta_from_zeta( this.zeta_from_z(z) );
        end
%% theta from zs ( zs --> z --> zeta --> [theta,r] )
        function [theta,r] = theta_from_zs(this,zs)
            [theta,r] = this.theta_from_z( this.scale_from_unit_chord(zs) );
        end

%% Complex Potential
        % function val = F_cylinder(this,zeta)
        %     uinf = 1;
        %     r    = this.a;
        %     gam  = 4*pi*uinf*r*sin( this.alpha + this.beta);
        %     eiam = exp(-1i*this.alpha);
        %     eiap = exp( 1i*this.alpha);
        %     z2  = zeta-this.mu;
        %     val = uinf*z2*eiam + 1i*(gam/(2*pi))*log(z2/(r*eiap)) + (uinf*(r^2)*eiap)./z2;
        % end
%% F from zeta ( zeta --> F )
        function F = F_cylinder(this,zeta)
            c1  = exp(-1i*this.alpha).*(zeta-this.mu);
            c2  = 2*1i*this.a*sin(this.alpha+this.beta)*log((zeta-this.mu)/ ( this.a*exp(1i*this.alpha) ) );
            c3  = this.a^2*exp(1i*this.alpha)./(zeta-this.mu);
            F = c1 + c2 + c3;
        end
%% psi from zeta ( zeta --> Im(F(zeta)) = psi )
        function psi = psi_from_zeta(this,zeta)
            psi = imag( this.F_cylinder( zeta ) );
        end
%% psi from z ( z --> zeta --> Im(F(zeta)) = psi )
        function psi = psi_from_z(this,z)
            psi = this.psi_from_zeta( this.zeta_from_z(z) );
        end
%% psi from zs ( zs --> z --> Im(F(zeta)) = psi )
        function psi = psi_from_zs(this,zs)
            psi = this.psi_from_z( this.scale_from_unit_chord(zs) );
        end
%% psi from [x,y] ( [x,y] --> zs --> z --> Im(F(zeta)) = psi )
        function psi = psi_from_xy(this,x,y)
            psi = this.psi_from_zs(x+1i*y);
        end

%% phi from zeta ( zeta --> Re(F(zeta)) = phi )
        function phi = phi_from_zeta(this,zeta)
            phi = real( this.F_cylinder( zeta ) );
        end
%% phi from z ( z --> zeta --> Re(F(zeta)) = phi )
        function phi = phi_from_z(this,z)
            phi = this.phi_from_zeta( this.zeta_from_z(z) );
        end
%% phi from zs ( zs --> z --> Re(F(zeta)) = phi )
        function phi = phi_from_zs(this,zs)
            phi = this.phi_from_z( this.scale_from_unit_chord(zs) );
        end
%% phi from [x,y] ( [x,y] --> zs --> z --> Re(F(zeta)) = phi )
        function phi = phi_from_xy(this,x,y)
            phi = this.phi_from_zs(x+1i*y);
        end
%% phi from [x,y] ( [x,y] --> zs --> z --> Re(F(zeta)) = phi ) (numeric)
        function phi = phi_from_xy_num(this,x,y)
            z = this.scale_from_unit_chord( x + 1i*y );
            phi = arrayfun( @(z) this.phi_from_zeta( this.zeta_from_z_num(z) ), z );
        end
%% y from [x,psi,y_guess] ( solve psi(x,y) == psi for y )
        function y = psi_y_from_x(this,psi,x,y_guess)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            obj_fun = @(x,y) this.psi_from_xy(x,y) - psi;
            y = arrayfun(@(x)fzero(@(y)obj_fun(x,y),y_guess,options),x);
        end
%% x from [y,phi,x_guess] ( solve phi(x,y) == phi for x )
        function x = phi_x_from_y(this,phi,y,x_guess)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            obj_fun = @(x,y) this.phi_from_xy(x,y) - phi;
            x = arrayfun(@(y)fzero(@(x)obj_fun(x,y),x_guess,options),y);
        end
%% [x,y] --> zs --> z --> zeta --> F[phi + 1i*psi]
        function F = phi_psi_from_xy(this,x,y)
            z = this.scale_from_unit_chord(x+1i*y);
            F = this.F_cylinder( this.zeta_from_z(z) );
        end
%% F --> zeta ( solve F(zeta) == F for zeta )
        function zeta = zeta_from_potential(this,F,zeta_guess)
            options = optimset('TolFun',1e-15,'TolX',1e-17,'Display','none');
            zeta = arrayfun( @(F) fsolve( @(zeta) F-this.F_cylinder( zeta ),...
                                            zeta_guess, options ), F );
        end
%% [phi,psi] --> F --> zeta --> z --> zs
        function z = z_from_phi_psi(this,phi,psi)
            F = phi + 1i*psi;
            options = optimset('TolFun',1e-15,'TolX',1e-17,'Display','none');
            z = arrayfun( @(F) fsolve( @(z) F-this.F_cylinder( ...
                                     this.zeta_from_z_num(z) ), ...
                                     ( 2*(imag(F)>=0) - 1 )*1i, options ), F );
            z = this.scale_to_unit_chord(z);
        end
%% [psi,theta] --> r
        function r = r_from_psi_theta(this,psi,theta)
            z1 = this.z_from_zeta( this.zeta_from_theta(theta));
            z1 = (z1-this.zLE)/this.chord;
            phi = this.phi_from_xy(real(z1),imag(z1));
            z2 = this.z_from_phi_psi(phi,psi);
            [~,r] = this.theta_from_zeta(this.zeta_from_z(this.chord*z2+this.zLE));
        end
%% minimum distance from airfoil
        function [d,theta] = distance_from_zs(airfoil,zs)
            theta0 = airfoil.theta_from_zs(zs);
            options = optimset('TolFun',1e-15,'TolX',1e-17,'Display','none');
            dist = @(theta) abs(zs-airfoil.zs_from_theta(theta));
            % theta = arrayfun( @(z) fminbnd( @(theta) dist(theta),1e-12,2*pi), zs );
            theta = arrayfun( @(z) fminunc( @(theta) dist(theta),theta0,options), zs );
            d = dist(theta);
        end

%% minimum distance between airfoil and streamline
        function [d,zs_airfoil,zs_streamline,theta,x] = airfoil_streamline_distance(this,psi)
            [theta,x,d] = arrayfun( @(psi) obj_fun_solve(this,psi), psi );
            zs_airfoil = this.zs_from_theta(theta);
            zs_streamline = x + 1i*this.psi_y_from_x(psi,x,d.*(2*(psi>=0)-1));
            function [theta,x,d] = obj_fun_solve(this,psi)
                x_guess = 0;
                y_guess = this.psi_y_from_x(psi,x_guess,2*(psi>=0)-1);
                theta_guess = this.theta_from_zs(0+1i*y_guess);
                f1 = @(x,psi) this.psi_y_from_x(psi,x,y_guess);
                f2 = @(theta) this.zs_from_theta(theta);
                objfun = @(x,psi) abs(f1(x(1),psi) - f2(x(2))).^2;
                options = optimset('TolFun',1e-15,'TolX',1e-17,'Display','none');
                x_opt = fminunc( @(x) objfun(x,psi),[x_guess,theta_guess], options );
                theta = x_opt(2);
                x     = x_opt(1);
                d     = objfun(x_opt,psi);
            end

        end


%% airfoil coordinates
        function [x,y] = output_airfoil_coords_theta(this,theta)
            z = this.zs_from_theta( theta );
            x = real(z);
            y = imag(z);
        end
%% airfoil coordinates ( parameterized w.r.t. function F ) (type 1)
        function [x,y,theta] = output_airfoil_coords1(this,N,F)
            t = this.arc_length_param2(linspace(0,1,N),F);
            % t = reparam_fun(1,@(t)this.z_from_zeta(this.zeta_from_theta(2*pi*t)),F,N,0,1);
            theta = 2*pi*t;
            [x,y] = output_airfoil_coords_theta(this,theta);
        end
%% airfoil coordinates ( parameterized w.r.t. function F ) (type 2)
        function [x,y,theta] = output_airfoil_coords2(this,N,F)
            t = linspace(0,1,N);
            theta = 2*pi*this.arc_length_param(F(t));
            [x,y] = output_airfoil_coords_theta(this,theta);
        end
        
%% Mapping Derivatives

%% d (zeta) / d (theta)
        function dzeta_dtheta = diff_zeta_from_theta(this,theta)
            % derivative of cylinder mapping w.r.t. theta
            dzeta_dtheta = 1i*this.a*exp(1i*(theta-this.beta));
        end
%% d^2 (zeta) / d (theta)^2
        function d2zeta_dtheta2 = diff2_zeta_from_theta(this,theta)
            % 2nd derivative of cylinder mapping w.r.t. theta
            d2zeta_dtheta2 = -this.a*exp(1i*(theta-this.beta));
        end
%% d (z) / d (zeta)
        function dz_dzeta = diff_z_from_zeta(this,zeta)
            % derivative of the transform w.r.t. zeta ( d/dzeta [z(zeta] )
            zeta_frac = ( (zeta - this.l)./(zeta + this.l) ).^this.n;
            factor = 4*(this.n*this.l)^2;
            dz_dzeta = factor*zeta_frac./( (zeta.^2-1).*(1-zeta_frac).^2 );
            dz_dzeta(isnan(dz_dzeta)) = 0;
            dz_dzeta(isinf(dz_dzeta)) = 0;
        end
%% d^2 (z) / d (zeta)^2
        function d2z_dzeta2 = diff2_z_from_zeta(this,zeta)
            zeta_p = (zeta+this.l).^this.n;
            zeta_m = (zeta-this.l).^this.n;
            factor = 8*(this.n*this.l)^2;
            num = factor*zeta_p.*zeta_m.* ...
                (zeta.*zeta_m - zeta.*zeta_p + this.l*this.n*(zeta_p + zeta_m) );
            den = (zeta_p - zeta_m).^3 .* (this.l^2 - zeta.^2).^2;
            d2z_dzeta2 = num./den;
            d2z_dzeta2(isnan(d2z_dzeta2)) = 0;
            d2z_dzeta2(isinf(d2z_dzeta2)) = 0;
        end

%% Inverse Mapping Derivatives

%% d (zeta) / d(z)
        function dzeta_dz = diff_zeta_from_z(this,z)
            % derivative of the inverse transform ( d/dz [zeta(z)] )
            %
            %                      2 /      2 l n  \1/n
            %                   4 l  | 1 - ------- |
            %  d zeta                \     z + l n /
            %  ------ = ------------------------------------------
            %    d z      2    2  2  / /    2 l n      \1/n     \2
            %           (z  - l  n ) | | - ------- + 1 |    - 1 |
            %                        \ \   z + l n     /        /
            ln = this.n*this.l;
            z_p = (z+ln);
            xn = 1/this.n;
            dzeta_dz = -( 4*this.l^2*(1 - (2*ln)./z_p).^xn )./ ...
                        ( (ln^2 - z.^2).*( (- (2*ln)./z_p + 1).^xn - 1 ).^2 );
        end

%% d^2 (zeta) / d (z)^2
        function d2zeta_dz2 = diff2_zeta_from_z(this,z)
            % 2nd derivative of the inverse transform ( d^2/dz^2 [zeta(z)] )
            % d2zeta_dz2 = -(this.l.^3.*(8.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) + 8.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n)) - 8.*this.l.^2.*z.*(((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) - ((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n)))./((this.l.^2.*this.n.^2 - z.^2).^2.*(3.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) - 3.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n) + ((z - this.l.*this.n)./(z + this.l.*this.n)).^(3./this.n) - 1));
            d2zeta_dz2 = -( this.l.^3.*( ...
                               8.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) ...
                             + 8.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n) ) ...
                             - 8.*this.l.^2.*z.*(((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) ...
                             - ((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n)) ) ...
                             ./ ( (this.l.^2.*this.n.^2 - z.^2).^2 ...
                                .*( 3.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) ...
                                  - 3.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n) ...
                                  + ((z - this.l.*this.n)./(z + this.l.*this.n)).^(3./this.n) - 1));
        end
%% Geometry
        function N = unit_normal_cmplx(this,theta)
            N = 1i*this.diff_z_from_zeta( this.zeta_from_theta(theta) ) ...
                     .*this.diff_zeta_from_theta(theta);
            N = -N; % outward pointing normal
            mag = max(abs(N),eps(1));
            N = N./mag;
        end
        function [N1,N2,N3] = unit_normal(this,theta)
            % outward pointing normal
            N = this.unit_normal_cmplx(theta);
            N1 = real(N);
            N2 = imag(N);
            N3 = 0*real(N);
        end
        function dS = airfoil_differential_arc_length(this,theta)
            dS = this.diff_z_from_zeta( this.zeta_from_theta(theta) ) ...
                     .*this.diff_zeta_from_theta(theta);
        end
        function K = airfoil_curvature1(this,theta)
            dz = this.diff_z_from_zeta( this.zeta_from_theta(theta) ) ...
                     .*this.diff_zeta_from_theta(theta);
            xdot = real(dz);
            ydot = imag(dz);
            % D[D[f1[f2[x]], x], x] = f2'[x]^2 f1''[f2[x]] + f1'[f2[x]] f2''[x]
            % dz2  = 
            term1 = this.diff_zeta_from_theta(theta).^2;
            term2 = this.diff2_z_from_zeta( this.zeta_from_theta(theta) );
            term3 = this.diff_z_from_zeta( this.zeta_from_theta(theta) );
            term4 = this.diff2_zeta_from_theta(theta);
            dF2 = term1.*term2 + term3.*term4;
            xdotdot = real(dF2);
            ydotdot = imag(dF2);
            K = abs(xdot.*ydotdot - ydot.*xdotdot)./(xdot.^2+ydot.^2).^(3/2);
        end
        function K = airfoil_curvature2(this,theta)
            dz = this.diff_z_from_zeta( this.zeta_from_theta(theta) ) ...
                     .*this.diff_zeta_from_theta(theta);
            % D[D[f1[f2[x]], x], x] = f2'[x]^2 f1''[f2[x]] + f1'[f2[x]] f2''[x]
            % dz2  = 
            term1 = this.diff_zeta_from_theta(theta).^2;
            term2 = this.diff2_z_from_zeta( this.zeta_from_theta(theta) );
            term3 = this.diff_z_from_zeta( this.zeta_from_theta(theta) );
            term4 = this.diff2_zeta_from_theta(theta);
            dz2 = term1.*term2 + term3.*term4;
            rdot_rdot       = real(dz.*conj(dz));
            rdot_rdotdot    = real(dz.*conj(dz2));
            rdotdot_rdotdot = real(dz2.*conj(dz2));
            K = sqrt(rdot_rdot.*rdotdot_rdotdot - rdot_rdotdot.^2)./(rdot_rdot).^(3/2);
        end
        function S = airfoil_arc_length(this,t0,t1)
            S = integral(@(t)diff_arc_length(this,t),t0,t1,"AbsTol",1e-15,'RelTol',1e-14);
            function dS = diff_arc_length(this,t)
                dS = 2*pi*abs(this.airfoil_differential_arc_length( 2*pi*t ));
                dS(isnan(dS)) = 0;
                dS(isinf(dS)) = 0;
            end
        end
        function ts = arc_length_param(this,t)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            N = numel(t);
            L_total = this.airfoil_arc_length(t(1),t(N));
            dL = ( t(N) - t(1) ) / L_total;
            ts = zeros(N,1);
            ts(1) = t(1);
            ts(N) = t(N);
            for i = 2:N-1
                ftmp = @(tsi) t(i)-t(i-1) - dL * this.airfoil_arc_length(ts(i-1),tsi);
                ts(i) = fzero( @(tsi)ftmp(tsi),ts(i-1),options);
            end
        end
        function ts = arc_length_param2(this,t,F)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            N = numel(t);
            L_total = this.airfoil_arc_length(t(1),t(N));
            dL = ( t(N) - t(1) ) / L_total;
            ts = zeros(N,1);
            ts(1) = t(1);
            ts(N) = t(N);
            for i = 2:N-1
                ftmp = @(tsi) F(t(i))-F(t(i-1)) - dL * this.airfoil_arc_length(ts(i-1),tsi);
                ts(i) = fzero( @(tsi)ftmp(tsi),ts(i-1),options);
            end
        end
        function bool = on_airfoil(this,x,y,scale,tol)
            z = x + 1i*y;
            if (scale)
                % z = z.*this.chord + this.airfoil_coords(this.thetaLE);
                z = this.scale_from_unit_chord(z);
            end
            zeta = this.zeta_from_z(z);
            zeta_dist = abs( zeta - this.mu );

            err = abs(zeta_dist - this.a);

            bool = err < tol;
        end
        function [t0,t1] = get_theta_from_z(this,x1,y1,x2,y2,scale)
            z1 = x1 + 1i*y1;
            z2 = x2 + 1i*y2;
            if (scale)
                % z1 = z1.*this.chord + this.airfoil_coords(this.thetaLE);
                % z2 = z2.*this.chord + this.airfoil_coords(this.thetaLE);
                z1 = this.scale_from_unit_chord(z1);
                z2 = this.scale_from_unit_chord(z2);
            end
            zeta1 = this.zeta_from_z(z1);
            zeta2 = this.zeta_from_z(z2);
            z_fun1 = @(theta) abs( this.zeta_from_theta(theta) - zeta1 );
            z_fun2 = @(theta) abs( this.zeta_from_theta(theta) - zeta2 );
            options = optimset('TolFun',1e-12,'TolX',1e-16);

            % single guess to avoid hopping the branch cut
            theta_guess = real( -1i*log( (zeta1-this.mu)/this.a ) );

            t0 = fminsearch(@(theta)z_fun1(theta),theta_guess,options);
            t1 = fminsearch(@(theta)z_fun2(theta),theta_guess,options);
        end
        function [z,dzdt,n] = curv_param(this,x1,y1,x2,y2,scale)
            [t0,t1] = this.get_theta_from_z(x1,y1,x2,y2,scale);
            map_fun = @(t) (t1-t0)*t + t0;
            dmap_fun_dt = @(t) t*0 + (t1-t0);
            z = @(t) this.z_from_theta(map_fun(t));
            
            dzdt = @(t) this.airfoil_differential_arc_length(map_fun(t)).*dmap_fun_dt(t);
            n = @(t) this.unit_normal_cmplx(map_fun(t));
        end
        function [z,dzdt,n] = line_param(this,x1,y1,x2,y2,scale)
            z1 = x1 + 1i*y1;
            z2 = x2 + 1i*y2;
            if (scale)
                % z1 = z1.*this.chord + this.airfoil_coords(this.thetaLE);
                % z2 = z2.*this.chord + this.airfoil_coords(this.thetaLE);
                z1 = this.scale_from_unit_chord(z1);
                z2 = this.scale_from_unit_chord(z2);
            end
            dzdt = @(t) t*0 + (z2-z1);
            z = @(t) (z2-z1).*t + z1;
            
            ntmp = -1i*(z2-z1);
            n = @(t) 0*t + ntmp./abs(ntmp);
        end
        function area = integrate_polygon_area(this,x,y,scale,use_curv)
            N = length(x);
            tol = 1e-8;
            area = 0;
            for i = 1:N
                ip1 = mod(i,N)+1;
                on_surface = ( this.on_airfoil( x(i),   y(i),   scale, tol ) ...;
                            && this.on_airfoil( x(ip1), y(ip1), scale, tol ) ...
                            && use_curv );
                if on_surface
                    [z,dzdt,~] = this.curv_param(x(i),y(i),x(ip1),y(ip1),scale);
                else
                    [z,dzdt,~] = this.line_param(x(i),y(i),x(ip1),y(ip1),scale);
                end
                zfun = @(t) conj(z(t)).*dzdt(t);
                area = area + integral(@(t)zfun(t),0,1,"AbsTol",1e-16,"RelTol",1e-12);
            end
            area = real( area/(2i) );
            if (scale)
                area = area/this.chord^2;
            end
        end
        function val = airfoil_surface_integral(this,fun,t0,t1)
            val = integral(@(theta)fun(theta) ...
                           .*abs(this.airfoil_differential_arc_length(theta)),...
                           t0,t1,"AbsTol",1e-14,'RelTol',1e-12);
        end

%% velocity around cylinder
        function val = w_cylinder(this,zeta)
            c1  = exp(-1i*this.alpha);
            c2  = 2*1i*this.a*sin(this.alpha+this.beta);
            c3  = -this.a^2*exp(1i*this.alpha);
            den = (zeta-this.mu);
            val = c1 + c2./den + c3./den.^2;
        end
        function val = dw_cylinder_dzeta(this,zeta)
            c2  = 2*1i*this.a*sin(this.alpha+this.beta);
            c3  = -this.a^2*exp(1i*this.alpha);
            den = (zeta-this.mu);
            val = -c2./den.^2 - 2*c3./den.^3;
        end
%% velocity around airfoil
        function val = w_airfoil(this,zeta)
            dz_dzeta = this.diff_z_from_zeta(zeta);
            wtilde   = this.w_cylinder(zeta);
            val      = this.vinf*wtilde./dz_dzeta;
        end
%% complex derivatives of velocity

%% d (w) / d (zeta)
        function val = dw_airfoil_dzeta(this,zeta)
            wtilde     = this.w_cylinder(zeta);
            dwtilde    = this.dw_cylinder_dzeta(zeta);
            dz_dzeta   = this.diff_z_from_zeta(zeta);
            d2z_dzeta2 = this.diff2_z_from_zeta(zeta);
            val        = this.vinf*( dwtilde./dz_dzeta ...
                                   - (wtilde./dz_dzeta.^2).*d2z_dzeta2 );
        end

%% d (w) / d (theta)
        function val = dw_airfoil_dtheta(this,theta)
            zeta         = this.zeta_from_theta(theta);
            dw_dzeta     = this.dw_airfoil_dzeta(this,zeta);
            dzeta_dtheta = this.diff_zeta_from_theta(theta);
            val          = dw_dzeta.*dzeta_dtheta;
        end

%% d (w) / d (z)
        function val = dw_airfoil_dz(this,z)
            zeta     = this.zeta_from_z(z);
            dzeta_dz = this.diff_zeta_from_z(z);
            dw_dzeta = this.dw_airfoil_dzeta(zeta);
            val      = dw_dzeta.*dzeta_dz;
        end

%% d (p) / d (z)
        function val = dp_airfoil_dz(this,z)
            zeta     = this.zeta_from_z(z);
            dzeta_dz = this.diff_zeta_from_z(z);
            dw_dzeta = this.dw_airfoil_dzeta(zeta);
            dp_dw    = -2*(1/2)*this.rhoinf*conj(this.w_airfoil(zeta));
            val      = dp_dw.*dw_dzeta.*dzeta_dz;
        end
%% d (p) / d (theta) (on surface of airfoil)
        function val = dp_airfoil_dtheta(this,theta)
            zeta         = this.zeta_from_theta(theta);
            dzeta_dtheta = this.diff_zeta_from_theta(theta);
            dw_dzeta     = this.dw_airfoil_dzeta(zeta);
            dp_dw        = -2*(1/2)*this.rhoinf*conj(this.w_airfoil(zeta));
            val          = dp_dw.*dw_dzeta.*dzeta_dtheta;
        end
%% Primitive Variables (function of zeta)
        function rho = airfoil_density(this,zeta)
            rho = this.rhoinf*ones(size(zeta));
        end
        function u = airfoil_x_velocity(this,zeta)
            u = real(this.w_airfoil(zeta));
        end
        function v = airfoil_y_velocity(this,zeta)
            v = -imag(this.w_airfoil(zeta));
        end
        function w = airfoil_z_velocity(~,zeta)
            w = zeros(size(zeta));
        end
        function vn = airfoil_contravariant_velocity(this,zeta,N1,N2,~)
            w  = this.w_airfoil(zeta);
            vn = real(w).*N1 - imag(w).*N2;
        end
        function p = airfoil_pressure(this,zeta)
            w2 = abs(this.w_airfoil(zeta)).^2;
            p  = this.pinf + (1/2)*this.rhoinf*( this.vinf^2 - w2 );
        end
%% Primitive Variables (function of [x,y])
        function V = get_prim_variables(this,i,x,y,scale)
            z = x + 1i*y;
            if (scale)
                % z = z.*this.chord + this.airfoil_coords(this.thetaLE);
                z = this.scale_from_unit_chord(x+1i*y);
            end
            zeta = this.zeta_from_z(z);
            switch(i)
                case(1)
                    V = this.airfoil_density(zeta);
                case(2)
                    V = this.airfoil_x_velocity(zeta);
                case(3)
                    V = this.airfoil_y_velocity(zeta);
                case(4)
                    V = this.airfoil_z_velocity(zeta);
                case(5)
                    V = this.airfoil_pressure(zeta);
                otherwise
                    error('only primitive variables 1-5 are supported')
            end
        end
%% Space Derivatives
        function val = dxvel_dx(this,z)
            val = real(this.dw_airfoil_dz(z));
        end
        function val = dxvel_dy(this,z)
            val = -imag(this.dw_airfoil_dz(z));
        end
        function val = dxvel_dz(~,z)
            val = zeros(size(z));
        end
        function val = dyvel_dx(this,z)
            val = this.dxvel_dy(z);
        end
        function val = dyvel_dy(this,z)
            val = -this.dxvel_dx(z);
        end
        function val = dyvel_dz(~,z)
            val = zeros(size(z));
        end
        function val = dzvel_dx(~,z)
            val = zeros(size(z));
        end
        function val = dzvel_dy(~,z)
            val = zeros(size(z));
        end
        function val = dzvel_dz(~,z)
            val = zeros(size(z));
        end
        function val = dp_dx(this,z)
            val = real(this.dp_airfoil_dz(z));
        end
        function val = dp_dy(this,z)
            val = -imag(this.dp_airfoil_dz(z));
        end
        function val = dp_dz(~,z)
            val = zeros(size(z));
        end
        function val = drho_dx(~,z)
            val = zeros(size(z));
        end
        function val = drho_dy(~,z)
            val = zeros(size(z));
        end
        function val = drho_dz(~,z)
            val = zeros(size(z));
        end
        function dVdxi = get_prim_variable_derivatives(this,i,d,x,y,scale)
            z = x + 1i*y;
            factor = 1;
            if (scale)
                z = z.*this.chord + this.z_from_theta(this.thetaLE);
                factor = this.chord;
            end
            switch(d)
                case(1)
                    switch(i)
                        case(1)
                            dVdxi = factor*this.drho_dx(z);
                        case(2)
                            dVdxi = factor*this.dxvel_dx(z);
                        case(3)
                            dVdxi = factor*this.dyvel_dx(z);
                        case(4)
                            dVdxi = factor*this.dzvel_dx(z);
                        case(5)
                            dVdxi = factor*this.dp_dx(z);
                        otherwise
                            error('only primitive variables 1-5 are supported')
                    end
                case(2)
                    switch(i)
                        case(1)
                            dVdxi = factor*this.drho_dy(z);
                        case(2)
                            dVdxi = factor*this.dxvel_dy(z);
                        case(3)
                            dVdxi = factor*this.dyvel_dy(z);
                        case(4)
                            dVdxi = factor*this.dzvel_dy(z);
                        case(5)
                            dVdxi = factor*this.dp_dy(z);
                        otherwise
                            error('only primitive variables 1-5 are supported')
                    end
                case(3)
                    switch(i)
                        case(1)
                            dVdxi = factor*this.drho_dz(z);
                        case(2)
                            dVdxi = factor*this.dxvel_dz(z);
                        case(3)
                            dVdxi = factor*this.dyvel_dz(z);
                        case(4)
                            dVdxi = factor*this.dzvel_dz(z);
                        case(5)
                            dVdxi = factor*this.dp_dz(z);
                        otherwise
                            error('only primitive variables 1-5 are supported')
                    end
            end
        end
        %% Conserved Variable Flux
        function F = mass_flux(this,zeta,N1,N2,~)
            w = this.w_airfoil(zeta);
            vn = real(w).*N1 - imag(w).*N2;
            F = this.airfoil_density(zeta).*vn;
        end
        function F = xmtm_flux(this,zeta,N1,N2,~)
            w   = this.w_airfoil(zeta);
            vn  = real(w).*N1 - imag(w).*N2;
            u   = real(w);
            rho = this.airfoil_density(zeta);
            p   = this.airfoil_pressure(zeta);
            F   = rho.*u.*vn + N1.*p;
        end
        function F = ymtm_flux(this,zeta,N1,N2,~)
            w = this.w_airfoil(zeta);
            vn = real(w).*N1 - imag(w).*N2;
            v   = -imag(w);
            rho = this.airfoil_density(zeta);
            p   = this.airfoil_pressure(zeta);
            F   = rho.*v.*vn + N2.*p;
        end
        function F = zmtm_flux(~,zeta,~,~,~)
            F = zeros(size(zeta));
        end
        function F = enrg_flux(this,zeta,N1,N2,~)
            w     = this.w_airfoil(zeta);
            vn    = real(w).*N1 - imag(w).*N2;
            v2    = abs(w).^2;
            p     = this.airfoil_pressure(zeta);
            rho   = this.airfoil_density(zeta);
            gxgm1 = this.gamma/(this.gamma-1);
            ht    = gxgm1*p./rho + 0.5*v2;
            F     = rho.*ht.*vn;
        end
        function F = mass_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.mass_flux(this.zeta_from_theta(theta),N1,N2,N3);
        end
        function F = xmtm_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.xmtm_flux(this.zeta_from_theta(theta),N1,N2,N3);
        end
        function F = ymtm_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.ymtm_flux(this.zeta_from_theta(theta),N1,N2,N3);
        end
        function F = zmtm_flux_on_airfoil(~,theta)
            F = zeros(size(theta));
        end
        function F = enrg_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.enrg_flux(this.zeta_from_theta(theta),N1,N2,N3);
        end
        function vals = integrate_fluxes_on_airfoil(this,t0,t1)
            fun1 = @(theta) this.mass_flux_on_airfoil(theta);
            fun2 = @(theta) this.xmtm_flux_on_airfoil(theta);
            fun3 = @(theta) this.ymtm_flux_on_airfoil(theta);
            fun5 = @(theta) this.enrg_flux_on_airfoil(theta);
            vals(1,1) = this.airfoil_surface_integral(fun1,t0,t1);
            vals(2,1) = this.airfoil_surface_integral(fun2,t0,t1);
            vals(3,1) = this.airfoil_surface_integral(fun3,t0,t1);
            vals(4,1) = 0;
            vals(5,1) = this.airfoil_surface_integral(fun5,t0,t1);
        end
        
        function vals = get_flux_integrals(this,x1,y1,x2,y2,scale,use_curv)
            tol = 1e-8;
            on_surface = this.on_airfoil(x1,y1,scale,tol) ...
                      && this.on_airfoil(x2,y2,scale,tol) ...
                      && use_curv;
            if on_surface
                [t0,t1] = this.get_theta_from_z(x1,y1,x2,y2,scale);
                vals = this.integrate_fluxes_on_airfoil(t0,t1);
            else
                [z,dzdt,normal] = this.line_param(x1,y1,x2,y2,scale);
                fun1 = @(t) this.mass_flux( this.zeta_from_z(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                fun2 = @(t) this.xmtm_flux( this.zeta_from_z(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                fun3 = @(t) this.ymtm_flux( this.zeta_from_z(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                fun5 = @(t) this.enrg_flux( this.zeta_from_z(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                vals(1,1) = integral(@(t)fun1(t),0,1,"AbsTol",1e-14,'RelTol',1e-12);
                vals(2,1) = integral(@(t)fun2(t),0,1,"AbsTol",1e-14,'RelTol',1e-12);
                vals(3,1) = integral(@(t)fun3(t),0,1,"AbsTol",1e-14,'RelTol',1e-12);
                vals(4,1) = 0;
                vals(5,1) = integral(@(t)fun5(t),0,1,"AbsTol",1e-14,'RelTol',1e-12);
            end
            if scale
                vals = vals/this.chord;
            end
        end
        %% Numerically integrated Force Coefficients
        function [Cx,Cy] = integrate_cp_on_airfoil(this,t0,t1)
            % integrates Cp from theta_0 to theta_1 on the surface, resolved
            % onto coordinate axes
            vals = this.integrate_fluxes_on_airfoil(t0,t1);
            Cx = -vals(2)./(0.5*this.rhoinf*this.vinf^2*this.chord);
            Cy = -vals(3)./(0.5*this.rhoinf*this.vinf^2*this.chord);
        end
        function [CL,CD] = get_CL_CD(this)
            % convert numerically integrate Cx,Cy to CL,CD
            [Cx,Cy] = this.integrate_cp_on_airfoil(0,2*pi);
            CL = Cy*cos(this.alpha) - Cx*sin(this.alpha);
            CD = Cy*sin(this.alpha) + Cx*cos(this.alpha);
        end
        function Cp = get_averaged_cp_on_segment(this,x1,y1,x2,y2,scale,use_curv)
            % compute average Cp on a segment specified by start and end
            % coordinates (used to get analytic results from a computational
            % grid)
            tol = 1e-8;
            point1_on_surface = this.on_airfoil(x1,y1,scale,tol);
            point2_on_surface = this.on_airfoil(x2,y2,scale,tol);
            if ~( point1_on_surface || point2_on_surface )
                fprintf('off surface\n');
            end
            on_surface = point1_on_surface ...
                      && point2_on_surface ...
                      && use_curv;
            if on_surface
                [t0,t1] = this.get_theta_from_z(x1,y1,x2,y2,scale);
                pfun = @(theta) ( this.airfoil_pressure(this.zeta_from_theta(theta)) - this.pinf )./(0.5*this.rhoinf*this.vinf^2);
                Cp_integral = airfoil_surface_integral(this,@(theta)pfun(theta),t0,t1)*sign(t1-t0);
                area = abs(airfoil_surface_integral(this,@(theta)0*theta+1,t0,t1));
                Cp = Cp_integral/area;
            else
                [z,dzdt,~] = this.line_param(x1,y1,x2,y2,scale);
                pfun = @(t) ( this.airfoil_pressure(this.zeta_from_z(z(t))) - this.pinf )./(0.5*this.rhoinf*this.vinf^2);
                dS   = @(t)abs(this.diff_z_from_zeta(z(t)).*dzdt(t));
                Cp_integral = integral(@(t)pfun(t).*dS(t),0,1,"AbsTol",1e-14,'RelTol',1e-12);
                % area = sqrt( (x2-x1).^2 + (y2-y1).^2);
                area = integral(@(t)dS(t),0,1,"AbsTol",1e-14,'RelTol',1e-12);
                Cp = Cp_integral/area;
            end
        end
        function Cp = get_averaged_cp_on_airfoil(this,x,y,scale,use_curv)
            N = length(x);
            Cp = zeros(N-1,1);
            for i = 1:N-1
                Cp(i) = get_averaged_cp_on_segment(this,x(i),y(i),x(i+1),y(i+1),scale,use_curv);
            end
        end
%% MMS source term
        function src = calculate_source(this,x,y,scale,use_curv)
            N = length(x);
            src = zeros(5,1);
            for i = 1:N
                ip1 = mod(i,N)+1;
                flux = this.get_flux_integrals(x(i),y(i),x(ip1),y(ip1),scale,use_curv);
                src = src + flux;
            end
            vol = this.integrate_polygon_area(x,y,scale,use_curv);
            src = src/vol;
        end
    end
end