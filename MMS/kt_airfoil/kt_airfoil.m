classdef kt_airfoil
    properties
        l(1,1)       double  = 1.0
        epsilon(1,1) double  = 0.0
        kappa(1,1)   double  = 0.0
        tau(1,1)     double  = 0.0
        rhoinf(1,1)  double  = 1.0
        vinf(1,1)    double  = 1.0
        pinf(1,1)    double  = 1.0
        gamma(1,1)   double  = 1.4
        chord(1,1)   double  = 1.0
        thetaLE(1,1) double  = 0.0
        thetaTE(1,1) double  = 0.0
        xLE(1,1)     double  = 0.0
        yLE(1,1)     double  = 0.0
        xTE(1,1)     double  = 0.0
        yTE(1,1)     double  = 0.0

        n(1,1)       double  = 2.0
        alpha(1,1)   double  = 0.0
        beta(1,1)    double  = 0.0
        a(1,1)       double  = 1.0
        mu(1,1)      double  = 0.0 + 1i*0.0
    end
    methods
        function this = kt_airfoil(epsilon,kappa,tau)
            this.epsilon = epsilon;
            this.kappa   = kappa;
            this.tau     = tau;
            this.n    = 2 - this.tau/pi;
            this.a    = this.l*sqrt( (1+this.epsilon)^2 + this.kappa^2 );
            this.mu   = this.l*(-this.epsilon + this.kappa*1i);
            this.beta = asin(this.l*this.kappa/this.a);
            this = this.get_airfoil_chord();
        end
        function V = get_prim_variables(this,i,x,y,scale)
            z = x + 1i*y;
            if (scale)
                z = z.*this.chord + this.airfoil_coords(this.thetaLE);
            end
            zeta = this.z_to_zeta(z);
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

        function V = get_prim_variable_derivatives(this,i,d,x,y,scale)
            z = x + 1i*y;
            factor = 1;
            if (scale)
                z = z.*this.chord + this.airfoil_coords(this.thetaLE);
                factor = this.chord;
            end
            switch(d)
                case(1)
                    switch(i)
                        case(1)
                            V = factor*this.ddens_dx(z);
                        case(2)
                            V = factor*this.dxvel_dx(z);
                        case(3)
                            V = factor*this.dyvel_dx(z);
                        case(4)
                            V = factor*this.dzvel_dx(z);
                        case(5)
                            V = factor*this.dpres_dx(z);
                        otherwise
                            error('only primitive variables 1-5 are supported')
                    end
                case(2)
                    switch(i)
                        case(1)
                            V = factor*this.ddens_dy(z);
                        case(2)
                            V = factor*this.dxvel_dy(z);
                        case(3)
                            V = factor*this.dyvel_dy(z);
                        case(4)
                            V = factor*this.dzvel_dy(z);
                        case(5)
                            V = factor*this.dpres_dy(z);
                        otherwise
                            error('only primitive variables 1-5 are supported')
                    end
                case(3)
                    switch(i)
                        case(1)
                            V = factor*this.ddens_dz(z);
                        case(2)
                            V = factor*this.dxvel_dz(z);
                        case(3)
                            V = factor*this.dyvel_dz(z);
                        case(4)
                            V = factor*this.dzvel_dz(z);
                        case(5)
                            V = factor*this.dpres_dz(z);
                        otherwise
                            error('only primitive variables 1-5 are supported')
                    end
            end
        end
        function this = set_alpha(this,alpha)
            this.alpha = deg2rad(alpha);
        end
        function val = CL(this)
            val = 8*pi*(this.a/this.chord)*sin(this.alpha+this.beta);
        end
        function val = CD(~)
            val = 0;
        end
        function val = CX(this)
            val = -this.CL * sin(this.alpha);
        end
        function val = CY(this)
            val = this.CL * cos(this.alpha);
        end

        function [h1,h2] = plot_airfoil(this,scale)
            if (scale)
                zfun = @(theta) ( this.airfoil_coords(theta) - this.airfoil_coords(this.thetaLE) )./this.chord;
            else
                zfun = @(theta) this.airfoil_coords(theta);
            end
            % h = fplot( @(theta)real(zfun(theta)), @(theta)imag(zfun(theta)), [0,2*pi] );
            h1 = fplot( @(theta)real(zfun(theta)), @(theta)imag(zfun(theta)), [0,this.thetaLE] );
            hold on;
            h2 = fplot( @(theta)real(zfun(theta)), @(theta)imag(zfun(theta)), [this.thetaLE,2*pi] );
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
                z = ( this.airfoil_coords(theta) - this.airfoil_coords(this.thetaLE) )./this.chord;
            else
                z = this.airfoil_coords(theta);
            end
            x = real(z(:)) - r(:).*N1(:);
            y = imag(z(:)) - r(:).*N2(:);
        end

        function [h1,h2] = plot_cp(this)
            tol = 1e-12;
            xfun = @(theta) real(this.airfoil_coords(theta) - this.airfoil_coords(this.thetaLE) )./this.chord;
            pfun = @(theta) ( this.airfoil_pressure(this.cylinder_map(theta)) - this.pinf )./(0.5*this.rhoinf*this.vinf^2);
            h1 = fplot( @(theta)xfun(theta), @(theta)pfun(theta), [0+tol,this.thetaLE] );
            hold on;
            h2 = fplot( @(theta)xfun(theta), @(theta)pfun(theta), [this.thetaLE,2*pi-tol] );
        end

        function zeta = cylinder_map(this,theta)
            zeta = this.a*exp(1i*(theta-this.beta)) + this.mu;
        end
        function dzeta_dtheta = cylinder_map_derivative(this,theta)
            dzeta_dtheta = 1i*this.a*exp(1i*(theta-this.beta));
        end
        function d2zeta_dtheta2 = cylinder_map_derivative2(this,theta)
            d2zeta_dtheta2 = -this.a*exp(1i*(theta-this.beta));
        end
        function z = zeta_to_z(this,zeta)
            zeta_p = (zeta+this.l).^this.n;
            zeta_m = (zeta-this.l).^this.n;
            z = this.n*this.l*( zeta_p + zeta_m )./( zeta_p - zeta_m );
        end
        function zeta = z_to_zeta(this,z)
            z_p = (z+this.n*this.l);
            z_m = (z-this.n*this.l);
            xn = 1/this.n;
            zeta = -this.l*( (z_m./z_p).^xn + 1) ./ ( (z_m./z_p).^xn - 1);
        end
        function dz_dzeta = diff_zeta_to_z(this,zeta)
            zeta_frac = ( (zeta - this.l)./(zeta + this.l) ).^this.n;
            factor = 4*(this.n*this.l)^2;
            dz_dzeta = factor*zeta_frac./( (zeta.^2-1).*(1-zeta_frac).^2 );
            dz_dzeta(isnan(dz_dzeta)) = 0;
            dz_dzeta(isinf(dz_dzeta)) = 0;
        end
        function dzeta_dz = diff_z_to_zeta(this,z)
            % -(4*this_l^2*(1 - (2*this_l*this_n)/(z + this_l*this_n))^(1/this_n))/((- z^2 + this_l^2*this_n^2)*((- (2*this_l*this_n)/(z + this_l*this_n) + 1)^(1/this_n) - 1)^2)
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
            dzeta_dz = -(4*this.l^2*(1 - (2*ln)./z_p).^xn)./((ln^2 - z.^2).*((- (2*ln)./z_p + 1).^xn - 1).^2);
            % dzeta_dz = (this.l.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n - 1).*((z - this.l.*this.n)./(z + this.l.*this.n).^2 - 1./(z + this.l.*this.n)))./(this.n.*(((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) - 1)) - (this.l.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n - 1).*((z - this.l.*this.n)./(z + this.l.*this.n).^2 - 1./(z + this.l.*this.n)).*(((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) + 1))./(this.n.*(((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) - 1).^2);
        end

        function d2zeta_dz2 = diff_z_to_zeta2(this,z)
            d2zeta_dz2 = -(this.l.^3.*(8.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) + 8.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n)) - 8.*this.l.^2.*z.*(((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) - ((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n)))./((this.l.^2.*this.n.^2 - z.^2).^2.*(3.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(1./this.n) - 3.*((z - this.l.*this.n)./(z + this.l.*this.n)).^(2./this.n) + ((z - this.l.*this.n)./(z + this.l.*this.n)).^(3./this.n) - 1));
        end

        function dz = diff_zeta_to_z2(this,zeta)
            % (8*l^2*n^2*(zeta + l)^n*(zeta - l)^n*(zeta*(zeta - l)^n - zeta*(zeta + l)^n + l*n*(zeta + l)^n + l*n*(zeta - l)^n))/(((l + zeta)^n - (- l + zeta)^n)^3*(l^2 - zeta^2)^2)
            % (8*l^2*n^2*zeta_p*zeta_m*(zeta*zeta_m - zeta*zeta_p + l*n*zeta_p + l*n*zeta_m))/((zeta_p - zeta_m)^3*(l^2 - zeta^2)^2)
            % (8*l^2*n^2*zeta_p*zeta_m*(zeta*zeta_m - zeta*zeta_p + l*n*zeta_p + l*n*zeta_m))
            zeta_p = (zeta+this.l).^this.n;
            zeta_m = (zeta-this.l).^this.n;
            factor = 8*(this.n*this.l)^2;
            num = factor*zeta_p.*zeta_m.*(zeta.*zeta_m - zeta.*zeta_p + this.l*this.n*(zeta_p + zeta_m) );
            den = (zeta_p - zeta_m).^3 .* (this.l^2 - zeta.^2).^2;
            dz = num./den;
            dz(isnan(dz)) = 0;
            dz(isinf(dz)) = 0;
        end
        
        function z = airfoil_coords(this,theta)
            z = this.zeta_to_z( this.cylinder_map(theta) );
        end
        function [x,y] = output_airfoil_coords_theta(this,theta)
            z = ( this.airfoil_coords(theta) - this.airfoil_coords(this.thetaLE) )./this.chord;
            x = real(z);
            y = imag(z);
        end
        function [x,y] = output_airfoil_coords1(this,N,F)
            t = this.arc_length_param(linspace(0,1,N));
            theta = 2*pi*F(t);
            [x,y] = output_airfoil_coords_theta(this,theta);
        end
        function [x,y] = output_airfoil_coords2(this,N,F)
            t = linspace(0,1,N);
            theta = 2*pi*this.arc_length_param(F(t));
            [x,y] = output_airfoil_coords_theta(this,theta);
        end
        function this = get_airfoil_chord(this)
            z_fun = @(s,theta)s*real(this.airfoil_coords(theta));
            this.thetaLE = fminbnd(@(theta)z_fun( 1,theta),0,2*pi);
            % this.thetaTE = fminbnd(@(theta)z_fun(-1,theta),0,2*pi);
            this.thetaTE = 0;
            zLE = this.airfoil_coords(this.thetaLE);
            zTE = this.airfoil_coords(this.thetaTE);
            this.chord = ( abs(zLE) + abs(zTE) );
            this.xLE   = real(zLE);
            this.yLE   = imag(zLE);
            this.xTE   = real(zTE);
            this.yTE   = imag(zTE);
        end
        function wtilde = cylinder_velocity(this,zeta)
            wtilde = exp(-1i*this.alpha) ...
                   + 2*1i*this.a*sin(this.alpha+this.beta)./(zeta-this.mu) ...
                   - this.a^2*exp(1i*this.alpha)./(zeta-this.mu).^2;
        end
        function dwtilde = cylinder_velocity_derivative(this,zeta)
            dwtilde = ...
                -2*1i*this.a*sin(this.alpha+this.beta)./(zeta-this.mu).^2 ...
                + 2*this.a^2*exp(1i*this.alpha)./(zeta-this.mu).^3;
        end
        function w = airfoil_velocity(this,zeta)
            dz = this.diff_zeta_to_z(zeta);
            w  = this.vinf*this.cylinder_velocity(zeta)./dz;
        end
        function dw = airfoil_velocity_derivative(this,zeta)
            wtilde  = this.cylinder_velocity(zeta);
            dwtilde = this.cylinder_velocity_derivative(zeta);
            dz      = this.diff_zeta_to_z(zeta);
            dz2     = this.diff_zeta_to_z2(zeta);
            dw = this.vinf*( dwtilde./dz - wtilde.*dz2./(dz.^2) );
        end
        function dw = airfoil_velocity_derivative_wrt_theta(this,theta)
            zeta    = this.cylinder_map(theta);
            dzeta   = this.cylinder_map_derivative(theta);
            wtilde  = this.cylinder_velocity(zeta);
            dwtilde = this.cylinder_velocity_derivative(zeta);
            dz      = this.diff_zeta_to_z(zeta);
            dz2     = this.diff_zeta_to_z2(zeta);
            dw = this.vinf*( ...
                             dwtilde.*dzeta./dz ...
                             - wtilde.*dz2.*dzeta./(dz.^2) );
        end
        function dw = dwdz(this,z)
            zeta     = this.z_to_zeta(z);
            dzeta_dz = this.diff_z_to_zeta(z);
            dw_dzeta = this.airfoil_velocity_derivative(zeta);
            dw = dzeta_dz.*dw_dzeta;
        end
        function dp = dpdz(this,z)
            zeta     = this.z_to_zeta(z);
            dzeta_dz = this.diff_z_to_zeta(z);
            dw_dzeta = this.airfoil_velocity_derivative(zeta);
            % dp_dw = -(1/2)*this.rhoinf*conj(this.airfoil_velocity(zeta));
            dp_dw = -this.rhoinf*conj(this.airfoil_velocity(zeta));
            dp = dzeta_dz.*dw_dzeta.*dp_dw;
        end
        % function dp = dpdz2(this,z)
        %     zeta     = this.z_to_zeta(z);
        %     dzeta_dz = this.diff_z_to_zeta(z);
        %     d2zeta_dz2 = this.diff_z_to_zeta2(z);
        %     dw_dzeta = this.airfoil_velocity_derivative(zeta);
        %     d2w_dzeta = this.airfoil_velocity_derivative2(zeta);
        %     dp_dw = -this.rhoinf*conj(this.airfoil_velocity(zeta));
        %     dp = dzeta_dz.*dw_dzeta.*dp_dw;
        % end
        %     % ux = vy, uy = -vx
        %     % f'(z) = ux + ivx = uy - ivy
        function d = dxvel_dx(this,z)
            d = real(this.dwdz(z));
        end
        function d = dxvel_dy(this,z)
            d = -imag(this.dwdz(z));
        end
        function d = dxvel_dz(~,z)
            d = zeros(size(z));
        end
        function d = dyvel_dx(this,z)
            d = this.dxvel_dy(z);
        end
        function d = dyvel_dy(this,z)
            d = -this.dxvel_dx(z);
        end
        function d = dyvel_dz(~,z)
            d = zeros(size(z));
        end
        function d = dzvel_dx(~,z)
            d = zeros(size(z));
        end
        function d = dzvel_dy(~,z)
            d = zeros(size(z));
        end
        function d = dzvel_dz(~,z)
            d = zeros(size(z));
        end

        function d = dpres_dx(this,z)
            d = real(this.dpdz(z));
        end
        function d = dpres_dy(this,z)
            d = -imag(this.dpdz(z));
        end
        function d = dpres_dz(~,z)
            d = zeros(size(z));
        end

        function d = ddens_dx(~,z)
            d = zeros(size(z));
        end
        function drho = ddens_dy(~,z)
            drho = zeros(size(z));
        end
        function drho = ddens_dz(~,z)
            drho = zeros(size(z));
        end

        function rho = airfoil_density(this,zeta)
            rho = this.rhoinf*ones(size(zeta));
        end
        function u = airfoil_x_velocity(this,zeta)
            u = real(this.airfoil_velocity(zeta));
        end
        function v = airfoil_y_velocity(this,zeta)
            v = -imag(this.airfoil_velocity(zeta));
        end
        function w = airfoil_z_velocity(~,zeta)
            w = zeros(size(zeta));
        end
        function vn = airfoil_contravariant_velocity(this,zeta,N1,N2,~)
            w = this.airfoil_velocity(zeta);
            vn = real(w).*N1 - imag(w).*N2;
        end
        function p = airfoil_pressure(this,zeta)
            p = this.pinf + (1/2)*this.rhoinf*( this.vinf^2 ...
                      - abs(this.airfoil_velocity(zeta)).^2 );
        end
        function F = mass_flux(this,zeta,N1,N2,~)
            w = this.airfoil_velocity(zeta);
            vn = real(w).*N1 - imag(w).*N2;
            F = this.airfoil_density(zeta).*vn;
        end
        function F = xmtm_flux(this,zeta,N1,N2,~)
            w = this.airfoil_velocity(zeta);
            vn = real(w).*N1 - imag(w).*N2;
            F = this.airfoil_density(zeta).*real(w).*vn ...
              + N1.*this.airfoil_pressure(zeta);
        end
        function F = ymtm_flux(this,zeta,N1,N2,~)
            w = this.airfoil_velocity(zeta);
            vn = real(w).*N1 - imag(w).*N2;
            F = this.airfoil_density(zeta).*(-imag(w)).*vn ...
              + N2.*this.airfoil_pressure(zeta);
        end
        function F = zmtm_flux(~,zeta,~,~,~)
            F = 0*real(zeta);
        end
        function F = enrg_flux(this,zeta,N1,N2,~)
            w = this.airfoil_velocity(zeta);
            vn = real(w).*N1 - imag(w).*N2;
            p   = this.airfoil_pressure(zeta);
            rho = this.rhoinf;
            gxgm1 = this.gamma/(this.gamma-1);
            ht = gxgm1*p/rho + 0.5*abs(w).^2;
            F = rho.*ht.*vn;
        end

        function F = mass_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.mass_flux(this.cylinder_map(theta),N1,N2,N3);
        end
        function F = xmtm_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.xmtm_flux(this.cylinder_map(theta),N1,N2,N3);
        end
        function F = ymtm_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.ymtm_flux(this.cylinder_map(theta),N1,N2,N3);
        end
        function F = zmtm_flux_on_airfoil(~,theta)
            F = 0*real(theta);
        end
        function F = enrg_flux_on_airfoil(this,theta)
            [N1,N2,N3] = this.unit_normal(theta);
            F = this.enrg_flux(this.cylinder_map(theta),N1,N2,N3);
        end
        function N = unit_normal_cmplx(this,theta)
            N = 1i*this.diff_zeta_to_z( this.cylinder_map(theta) ) ...
                     .*this.cylinder_map_derivative(theta);
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
            dS = this.diff_zeta_to_z( this.cylinder_map(theta) ) ...
                     .*this.cylinder_map_derivative(theta);
        end
        function K = airfoil_curvature1(this,theta)
            dz = this.diff_zeta_to_z( this.cylinder_map(theta) ) ...
                     .*this.cylinder_map_derivative(theta);
            xdot = real(dz);
            ydot = imag(dz);
            % D[D[f1[f2[x]], x], x] = f2'[x]^2 f1''[f2[x]] + f1'[f2[x]] f2''[x]
            % dz2  = 
            term1 = this.cylinder_map_derivative(theta).^2;
            term2 = this.diff_zeta_to_z2( this.cylinder_map(theta) );
            term3 = this.diff_zeta_to_z( this.cylinder_map(theta) );
            term4 = this.cylinder_map_derivative2(theta);
            dF2 = term1.*term2 + term3.*term4;
            xdotdot = real(dF2);
            ydotdot = imag(dF2);
            K = abs(xdot.*ydotdot - ydot.*xdotdot)./(xdot.^2+ydot.^2).^(3/2);
        end
        function K = airfoil_curvature2(this,theta)
            dz = this.diff_zeta_to_z( this.cylinder_map(theta) ) ...
                     .*this.cylinder_map_derivative(theta);
            % D[D[f1[f2[x]], x], x] = f2'[x]^2 f1''[f2[x]] + f1'[f2[x]] f2''[x]
            % dz2  = 
            term1 = this.cylinder_map_derivative(theta).^2;
            term2 = this.diff_zeta_to_z2( this.cylinder_map(theta) );
            term3 = this.diff_zeta_to_z( this.cylinder_map(theta) );
            term4 = this.cylinder_map_derivative2(theta);
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
        function val = airfoil_surface_integral(this,fun,t0,t1)
            val = integral(@(theta)fun(theta) ...
                           .*abs(this.airfoil_differential_arc_length(theta)),t0,t1,"AbsTol",1e-14,'RelTol',1e-12);
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
        function [Cx,Cy] = integrate_cp_on_airfoil(this,t0,t1)
            vals = this.integrate_fluxes_on_airfoil(t0,t1);
            Cx = -vals(2)./(0.5*this.rhoinf*this.vinf^2*this.chord);
            Cy = -vals(3)./(0.5*this.rhoinf*this.vinf^2*this.chord);
        end
        function [CL,CD] = get_CL_CD(this)
            [Cx,Cy] = this.integrate_cp_on_airfoil(0,2*pi);
            CL = Cy*cos(this.alpha) - Cx*sin(this.alpha);
            CD = Cy*sin(this.alpha) + Cx*cos(this.alpha);
        end
        function bool = on_airfoil(this,x,y,scale,tol)
            z = x + 1i*y;
            if (scale)
                z = z.*this.chord + this.airfoil_coords(this.thetaLE);
            end
            zeta = this.z_to_zeta(z);
            zeta_dist = abs( zeta - this.mu );

            err = abs(zeta_dist - this.a);

            bool = err < tol;
        end
        function [t0,t1] = get_theta_from_z(this,x1,y1,x2,y2,scale)
            z1 = x1 + 1i*y1;
            z2 = x2 + 1i*y2;
            if (scale)
                z1 = z1.*this.chord + this.airfoil_coords(this.thetaLE);
                z2 = z2.*this.chord + this.airfoil_coords(this.thetaLE);
            end
            zeta1 = this.z_to_zeta(z1);
            zeta2 = this.z_to_zeta(z2);
            z_fun1 = @(theta) abs( this.cylinder_map(theta) - zeta1 );
            z_fun2 = @(theta) abs( this.cylinder_map(theta) - zeta2 );
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
            z = @(t) this.airfoil_coords(map_fun(t));
            
            dzdt = @(t) this.airfoil_differential_arc_length(map_fun(t)).*dmap_fun_dt(t);
            n = @(t) this.unit_normal_cmplx(map_fun(t));
        end
        function [z,dzdt,n] = line_param(this,x1,y1,x2,y2,scale)
            z1 = x1 + 1i*y1;
            z2 = x2 + 1i*y2;
            if (scale)
                z1 = z1.*this.chord + this.airfoil_coords(this.thetaLE);
                z2 = z2.*this.chord + this.airfoil_coords(this.thetaLE);
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
                fun1 = @(t) this.mass_flux( this.z_to_zeta(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                fun2 = @(t) this.xmtm_flux( this.z_to_zeta(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                fun3 = @(t) this.ymtm_flux( this.z_to_zeta(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
                fun5 = @(t) this.enrg_flux( this.z_to_zeta(z(t) ),real(normal(t)),imag(normal(t)),0).*abs(dzdt(t));
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
        function Cp = get_averaged_cp_on_segment(this,x1,y1,x2,y2,scale,use_curv)
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
                pfun = @(theta) ( this.airfoil_pressure(this.cylinder_map(theta)) - this.pinf )./(0.5*this.rhoinf*this.vinf^2);
                Cp_integral = airfoil_surface_integral(this,@(theta)pfun(theta),t0,t1)*sign(t1-t0);
                area = abs(airfoil_surface_integral(this,@(theta)0*theta+1,t0,t1));
                Cp = Cp_integral/area;
            else
                [z,dzdt,~] = this.line_param(x1,y1,x2,y2,scale);
                pfun = @(t) ( this.airfoil_pressure(this.z_to_zeta(z(t))) - this.pinf )./(0.5*this.rhoinf*this.vinf^2);
                dS   = @(t)abs(this.diff_zeta_to_z(z(t)).*dzdt(t));
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

        
    end
end