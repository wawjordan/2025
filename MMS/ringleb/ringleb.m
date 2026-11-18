classdef ringleb
    properties
        q_min(1,1)       double = 68
        q_max(1,1)
        k_min(1,1)       double = 0.7
        k_max(1,1)       double = 1.5
        rho0(1,1)        double = 1.0
        p0(1,1)          double = 100000.0
        a0(1,1)
        gamma(1,1)       double = 1.4 % constant
        rho_ref(1,1)     double = 1.0
        a_ref(1,1)       double = 340.0
        l_ref(1,1)       double = 1.0
        l_ref_grid(1,1)  double = 1.0
    end
    methods
%% Constructor
        function this = ringleb(varargin)
            validScalarNum    = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
            p = inputParser;
            p.addOptional('q_min',68,validScalarPosNum);
            p.addOptional('k_min',0.6,validScalarPosNum);
            p.addOptional('k_max',1.25,validScalarPosNum);
            p.addOptional('rho0',1.0,validScalarPosNum);
            p.addOptional('p0',100000.0,validScalarPosNum);
            p.addOptional('rho_ref',1.0,validScalarPosNum)
            p.addOptional('a_ref',340.0,validScalarPosNum)
            p.addOptional('l_ref',1.0,validScalarPosNum)
            p.addOptional('l_ref_grid',1.0,validScalarPosNum)
            parse(p,varargin{:});
            this.l_ref       = p.Results.l_ref;
            this.l_ref_grid  = p.Results.l_ref_grid;
            % L = this.l_ref / this.l_ref_grid;
            this.rho_ref     = p.Results.rho_ref;
            this.a_ref       = p.Results.a_ref;
            this.rho0        = p.Results.rho0 / this.rho_ref;
            this.p0          = p.Results.p0 ...
                             / (this.rho_ref * this.a_ref^2);
            this.a0        = sqrt( this.gamma * this.p0 / this.rho0 );
            this.q_max     = (4/(this.gamma+1));
            this.q_min     = p.Results.q_min / this.a_ref;
            this.k_min     = p.Results.k_min;
            this.k_max     = p.Results.k_max;
        end
        function grid = make_grid(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            psi = linspace(1./this.k_max,1./this.k_min,N_j);
            theta0 = asin(psi(end)*this.q_min);
            phi_max = this.phi_from_q_theta(this.q_min,theta0);
            phi_min = 0;
            phi = linspace(phi_min,phi_max,N_i);
            [PHI,PSI] = ndgrid(phi,psi);
            [X,Y] = x_y_from_phi_psi(this,PHI,PSI);
            grid.x = X;
            grid.y = Y;
            function [x,y] = x_y_from_phi_psi(this,phi,psi)
                q0 = 1/psi(1);
                q = arrayfun(@(psi,phi)fzero( @(q)obj_fun(q,1./psi,phi),...
                                              q0 ), psi, phi );
                rho = this.rho_value(q);
                J = this.J_value(q);
                x = 1.0./rho .* (0.5./(q.^2) - psi.^2) + 0.5*J;
                y = psi./(rho.*q) .* sqrt( max(1.0 - (q.*psi).^2,0.0));
                function f = obj_fun(q,k,phi)
                    f = q.^2 .* (phi.^2 .*( 1- 0.2*q.^2).^5 + 1./k.^2) - 1;
                end
            end
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
                                        [this.q_min,this.q_max],opts),x,y);
            psi = this.psi_from_x_q(x,q);
            arg = psi.*q;
            arg(arg<-1) = -1;
            arg(arg>1)  =  1;
            theta = asin(arg);
            function f = obj_fun(this,x,y,q)
                J = this.J_value(q);
                r = this.rho_value(q);
                f = 2.0*r*sqrt((x-0.5*J).^2 + y.^2);
                f = sqrt(1.0./f);
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
        function p = pressure_value(this,q)
            c = this.c_value(q);
            p = c.^(2.0*this.gamma/(this.gamma-1.0));
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