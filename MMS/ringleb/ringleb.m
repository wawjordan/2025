classdef ringleb
    properties
        v_min(1,1)       double = 68
        psi_min(1,1)     double = 1.25
        psi_max(1,1)     double = 4.0
        rho0(1,1)        double = 1.0
        p0(1,1)          double = 100000.0
        a0(1,1)          double = 340.0
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
            p.addOptional('v_min',68,validScalarPosNum);
            p.addOptional('psi_min',1.25,validScalarPosNum);
            p.addOptional('psi_max',4.0,validScalarPosNum);
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
            this.v_min       = p.Results.v_min / this.a_ref;
            this.rho0        = p.Results.rho0 / this.rho_ref;
            this.p0          = p.Results.p0 ...
                             / (this.rho_ref * this.a_ref^2);
            this.psi_min     = p.Results.psi_min;
            this.psi_max     = p.Results.psi_max;
        end
        function grid = make_grid(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            psi_0 = linspace( this.psi_min, this.psi_max, N_j );
            V   = zeros(N_i,N_j);
            psi = zeros(N_i,N_j);

            for i = 1:N_i
                psi(i,:) = psi_0;
            end

            for j = 1:N_j
                theta_min = asin(this.v_min*psi_0(j));
                phi_start = 1./(psi_0(j)*this.v_min)*cos(theta_min);
                phi_end   = -phi_start;
                phi = linspace(phi_start,phi_end,N_i);
                for i = 1:ceil((N_i+1)/2)
                    k = 1/psi_0(j);
                    V(i,j) = sqrt(k^2/(1+phi(i)^2));
                end
            end
            grid.x = ringleb.x_value(psi,V);
            grid.y = ringleb.y_value(psi,V);
        end
        function grid = make_grid2(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            % v = sin(theta.')./this.psi_min;
            % [THETA,V] = ndgrid(theta,v);
            % PSI = sin(THETA)./V;
            psi = linspace(this.psi_min,this.psi_max,N_j);
            theta0 = asin(psi*this.v_min);
            t     = linspace( 0,1,N_i).';
            theta = theta0 + t.*(pi/2 - theta0);
            PSI = repmat(psi,N_i,1);
            V = sin(theta)./PSI;
            
            grid.x = ringleb.x_value(PSI,V);
            grid.y = ringleb.y_value(PSI,V);

            % grid1 = this.make_grid(N_i,N_j);
            % grid.x = 0.5*(grid.x + grid1.x);
            % grid.y = 0.5*(grid.y + grid1.y);
        end
        function grid = make_grid3(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            % psi = 1./linspace(1./this.psi_min,1./this.psi_max,N_j);
            psi = linspace(this.psi_min,this.psi_max,N_j);
            theta0 = asin(psi*this.v_min);
            phi_max = ringleb.phi_value(theta0(1),this.v_min);
            phi_min = 0;
            phi = linspace(phi_min,phi_max,N_i);
            [PHI,PSI] = ndgrid(phi,psi);
            Z = z_from_phi_psi(PHI,PSI);
            grid.x = real(Z);
            grid.y = imag(Z);
            % function z = z_from_q_k(q,k)
            %     b   = ringleb.b_value(q); % b = sqrt( 1 - (gamma-1)/2 * V^2 )
            %     rho = b.^5;               % rho = b^(2/(gamma-1))
            %     L   = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) ...
            %         - 0.5*log((1.0+b)./(1.0-b));
            %     x = 1.0./rho .* (0.5./(q.^2) - 1./k.^2) + 0.5*L;
            %     y = 1./(k.*rho.*q) .* sqrt( 1.0 - (q./k).^2);
            %     z = x + 1i*y;
            % end
            function z = z_from_phi_psi(phi,psi)
                q0 = 1/psi(1);
                q = arrayfun(@(psi,phi)fzero(@(q)obj_fun(q,1./psi,phi),q0),psi,phi);
                b   = ringleb.b_value(q); % b = sqrt( 1 - (gamma-1)/2 * V^2 )
                rho = b.^5;               % rho = b^(2/(gamma-1))
                L   = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) ...
                    - 0.5*log((1.0+b)./(1.0-b));
                x = 1.0./rho .* (0.5./(q.^2) - psi.^2) + 0.5*L;
                y = psi./(rho.*q) .* sqrt( 1.0 - (q.*psi).^2);
                z = x + 1i*y;
                function f = obj_fun(q,k,phi)
                    f = q.^2 .* ( phi.^2 .*( 1- 0.2*q.^2).^5 + 1./k.^2) - 1;
                end
            end
        end

        function rho = density(this,x,y)
            [rho,~,~,~] = this.primitive_vars(x,y);
        end
        function u = x_velocity(this,x,y)
            [~,u,~,~] = this.primitive_vars(x,y);
        end
        function v = y_velocity(this,x,y)
            [~,~,v,~] = this.primitive_vars(x,y);
        end
        function p = pressure(this,x,y)
            [~,~,~,p] = this.primitive_vars(x,y);
        end
        function [V,psi] = v_psi_from_x_y(this,x,y)
            V   = ringleb.V_value(x,y,this.v_min,this.psi_min);
            psi = ringleb.psi_value(x,V);
        end

        function Vpsi = v_psi_from_z(this,z)
            V   = ringleb.V_value(real(z),imag(z),this.v_min,this.psi_min);
            psi = ringleb.psi_value(real(z),V);
            Vpsi = complex(V,psi);
        end
        function [rho,u,v,p] = primitive_vars(this,x,y)
            V   = ringleb.V_value(x,y,this.v_min,this.psi_min);
            b   = ringleb.b_value(V); % b = sqrt( 1 - (gamma-1)/2 * V^2 )
            rho = b.^5;               % rho = b^(2/(gamma-1))
            p   = b.^7 / 1.4;         % p = rho*b^2 / gamma
            L   = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) ...
                - 0.5*log((1.0+b)./(1.0-b));
            psi   = sqrt(0.5./V.^2 - rho.*(x-0.5*L));
            arg = psi.*V;
            arg(arg<-1) = -1;
            arg(arg>1)  =  1;
            theta = asin(arg);
            u   = -V*cos(theta);
            v   = -V*sin(theta);

            rho = rho*this.rho0;
            p   = p*this.p0;
            u   = u*this.a0;
            v   = v*this.a0;
        end
    end
    methods (Static)
        function b = b_value(V)
            b = sqrt(1.0 - 0.2*V.^2);
        end
        function L = L_value(V)
            b = ringleb.b_value(V);
            L = 1.0./b + 1.0./(3.0*b.^3) + 1.0./(5.0*b.^5) ...
                       - 0.5*log((1.0+b)./(1.0-b)); 
        end
        function rho = density_value(V)
            b = ringleb.b_value(V);
            rho = b.^5;
        end
        function rho = pressure_value(V)
            b = ringleb.b_value(V);
            rho = b.^7;
        end
        function x = x_value(psi,V)
            rho = ringleb.density_value(V);
            L   = ringleb.L_value(V);
            x = 1.0./rho .* (0.5./(V.^2) - psi.^2) + 0.5*L;
        end
        function y = y_value(psi,V)
            rho = ringleb.density_value(V);
            rad = 1.0 - (V.^2).*(psi.^2);
            rad = (rad>=0).*rad + (rad<0)*0;
            y = psi./(rho.*V) .* sqrt(rad);
        end
        function psi = psi_value(x,V)
            rho = ringleb.density_value(V);
            L   = ringleb.L_value(V);
            psi = sqrt(0.5./(V.^2) - rho.*(x-0.5*L));
        end
        function phi = phi_value(theta,V)
            rho = ringleb.density_value(V);
            phi = cos(theta)./(V.*rho);
        end

        function k = k_from_V_phi(phi,V)
            rho = ringleb.density_value(V);
            k = V./sqrt(1-phi.^2.*V.^2*rho.^2);
        end
        function theta = theta_value(psi,V)
            arg = psi.*V;
            arg(arg<-1) = -1;
            arg(arg>1)  =  1;
            theta = asin(arg);
        end
        function V = V_fun(x,y,V)
            rho = ringleb.density_value(V);
            L   = ringleb.L_value(V);
            den = 2.0*rho*sqrt((x-0.5*L).^2 + y.^2);
            V = sqrt(1.0./den);
        end
        function V = V_value(x,y,v_min,psi_min)
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            V0 = 0.5*(v_min+1.0/psi_min);
            V = arrayfun(@(x,y)fzero(@(V_)ringleb.V_fun(x,y,V_),V0,options),x,y);
        end
    end
end