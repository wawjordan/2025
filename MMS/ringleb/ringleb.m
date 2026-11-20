classdef ringleb
    properties
        q_min(1,1)       double = 68  % min. velocity [m/s]
        q_max(1,1)                    % max. velocity [m/s]
        k_min(1,1)       double = 0.6 % streamline parameter [nondim]
        k_max(1,1)       double = 1.5 % streamline parameter [nondim]
        rho0(1,1)        double = 1.0 % stagnation density [kg/m^3]
        p0(1,1)          double = 100000.0 % stagnation pressure [Pa]
        a0(1,1)                       % speed of sound at stagnation [m/s]
        gamma(1,1)       double = 1.4 % needs to be 1.4 for closed form
        rho_ref(1,1)     double = 1.0 % reference density [kg/m^3]
        a_ref(1,1)       double = 340.0 % reference speed of sound [m/s]
        l_ref(1,1)       double = 1.0   % reference length [m]
        l_ref_grid(1,1)  double = 1.0   % grid scale [nondim]
    end
    methods
%% Constructor
        function this = ringleb(varargin)
            q_min_default      = 68;
            k_min_default      = 0.6;
            k_max_default      = 1.25;
            rho0_default       = 1.0;
            p0_default         = 100000.0;
            rho_ref_default    = 1.0;
            a_ref_default      = 340.0;
            l_ref_default      = 1.0;
            l_ref_grid_default = 1.0;
            validScalarNum    = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
            p = inputParser;
            p.addOptional('q_min', q_min_default, validScalarPosNum);
            p.addOptional('k_min', k_min_default, validScalarPosNum);
            p.addOptional('k_max', k_max_default, validScalarPosNum);
            p.addOptional('rho0',  rho0_default,  validScalarPosNum);
            p.addOptional('p0',    p0_default,    validScalarPosNum);
            p.addOptional('rho_ref', rho_ref_default, validScalarPosNum);
            p.addOptional('a_ref',   a_ref_default,   validScalarPosNum);
            p.addOptional('l_ref',   l_ref_default,   validScalarPosNum);
            p.addOptional('l_ref_grid', l_ref_grid_default, validScalarPosNum);
            parse(p,varargin{:});
            this.l_ref       = p.Results.l_ref;
            this.l_ref_grid  = p.Results.l_ref_grid;
            % L = this.l_ref / this.l_ref_grid;
            % reference quantities remain dimensional
            this.rho_ref     = p.Results.rho_ref;
            this.a_ref       = p.Results.a_ref;

            % non-dimensional maximum possible velocity (limiting line)
            this.q_max     = (4/(this.gamma+1));

            % already non-dimensional
            this.k_min     = p.Results.k_min;
            this.k_max     = p.Results.k_max;

            % non-dimensionalize inputs
            this.rho0        = p.Results.rho0 / this.rho_ref;
            this.p0          = p.Results.p0 ...
                             / (this.rho_ref * this.a_ref^2);
            this.a0        = sqrt( this.gamma * this.p0 / this.rho0 );
            this.q_min     = p.Results.q_min / this.a_ref;
            
        end
        function grid = make_grid_phi_psi(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            psi = linspace(1./this.k_max,1./this.k_min,N_j);
            theta0 = asin(psi(end)*this.q_min);
            phi_max = this.phi_from_q_theta(this.q_min,theta0);
            phi_min = 0;
            phi = linspace(phi_min,phi_max,N_i);
            [PHI,PSI] = ndgrid(phi,psi);
            [X,Y] = this.xy_from_phi_psi(PHI,PSI);
            grid.x = X;
            grid.y = Y;
        end
        function grid = make_grid_k_q(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            psi_0 = 1./linspace(this.k_max,this.k_min,N_j);
            Q = zeros(N_i,N_j);
            PSI = zeros(N_i,N_j);
            for i = 1:N_i
                PSI(i,:) = psi_0;
            end

            for j = 1:N_j
                theta_min = asin(psi_0(j)*this.q_min);
                theta_max = pi - asin(psi_0(j)*this.q_min);
                theta = linspace(theta_min,theta_max,N_i);
                Q(:,j) = sin(theta)/psi_0(j);
            end
            [X,Y] = this.xy_from_q_psi(Q,PSI);
            grid.x = X;
            grid.y = Y;
        end
        function [x,y] = xy_from_phi_psi(this,phi,psi)
            q0 = 1/psi(1);
            q = arrayfun(@(psi,phi)fzero( @(q)obj_fun(q,1./psi,phi),...
                q0 ), psi, phi );
            [x,y] = this.xy_from_q_psi(q,psi);
            function f = obj_fun(q,k,phi)
                f = q.^2 .* (phi.^2 .*( 1- 0.2*q.^2).^5 + 1./k.^2) - 1;
            end
        end
        function [x,y] = xy_from_q_psi(this,q,psi)
            rho = this.rho_value(q);
            J = this.J_value(q);
            x = 1.0./rho .* (0.5./(q.^2) - psi.^2) + 0.5*J;
            y = psi./(rho.*q) .* sqrt( max(1.0 - (q.*psi).^2,0.0));
        end

        function rho = density(this,x,y)
            [rho,~,~,~] = this.primitive_vars_from_x_y(x,y);
        end
        function u = x_velocity(this,x,y)
            [~,u,~,~] = this.primitive_vars_from_x_y(x,y);
        end
        function v = y_velocity(this,x,y)
            [~,~,v,~] = this.primitive_vars_from_x_y(x,y);
        end
        function p = pressure(this,x,y)
            [~,~,~,p] = this.primitive_vars_from_x_y(x,y);
        end
        function [rho,u,v,p] = primitive_vars_from_x_y(this,x,y)
            [q,theta] = this.q_theta_from_x_y(x,y);
            [rho,u,v,p] = this.primitive_vars_q_theta(q,theta);
        end
        function [q,theta] = q_theta_from_x_y(this,x,y)
            opts =  optimset('TolFun',1e-15,'TolX',1e-17);
            q = arrayfun( @(x,y) fzero( @(q)obj_fun(this,x,y,q),...
                                        [1e-6,this.q_max],opts),x,y);
            psi = this.psi_from_x_q(x,q);
            arg = psi.*q;
            arg(arg<-1) = -1;
            arg(arg>1)  =  1;
            theta = asin(arg);
            function f = obj_fun(this,x,y,q)
                J = this.J_value(q);
                r = this.rho_value(q);
                f = 2.0*r*sqrt((x-0.5*J).^2 + y.^2);
                f = q-sqrt(1.0./f);
            end
        end
        function [rho,u,v,p] = primitive_vars_q_theta(this,q,theta)
            rho = this.rho_value(q)*this.rho0;
            p   = this.pressure_value(q)*this.p0;
            u   = -q.*cos(theta)*this.a0;
            v   = -q.*sin(theta)*this.a0;
        end
        function c = c_value(this,q)
            c = sqrt(1.0 - ((this.gamma-1.0)/2.0)*q.^2);
        end
        function J = J_value(this,q)
            c = this.c_value(q);
            J = 1.0./c + 1.0./(3.0*c.^3) + 1.0./(5.0*c.^5) ...
              - 0.5*log((1.0+c)./(1.0-c));
        end
        function rho = rho_value(this,q)
            c = this.c_value(q);
            rho = c.^(2.0/(this.gamma-1.0));
        end
        function M = mach_number(this,q)
            p   = this.pressure_value(q);
            rho = this.rho_value(q);
            M = q./sqrt(this.gamma*p./rho);
        end
        function M = mach_number2(this,q)
            M = sqrt(q.^2./(1-0.5*(this.gamma-1)*q.^2));
        end
        function p = pressure_value(this,q)
            c = this.c_value(q);
            p = c.^(2.0*this.gamma/(this.gamma-1.0))/this.gamma;
        end
        function psi = psi_from_x_q(this,x,q)
            rho = this.rho_value(q);
            J   = this.J_value(q);
            psi = sqrt(0.5./(q.^2) - rho.*(x-0.5*J));
        end
        function phi = phi_from_q_theta(this,q,theta)
            rho = this.rho_value(q);
            phi = cos(theta)./(q.*rho);
        end
    end
end