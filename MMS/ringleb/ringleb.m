classdef ringleb
    properties
        q_min(1,1)       double = 68  % min. velocity [m/s]
        q_max(1,1)                    % max. velocity [m/s]
        k_min(1,1)       double = 0.6 % streamline parameter [nondim]
        k_max(1,1)       double = 1.5 % streamline parameter [nondim]
        xk_min(1,1)      double =-1.0 % location on x-axis of k_min [m]
        xk_max(1,1)      double = 0.0 % location on x-axis of k_max [m]
        xq_min(1,1)      double = 1.0 % x-location of point for q_min,k_max
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
            xk_min_default     = NaN;
            xk_max_default     = NaN;
            xq_min_default     = NaN;
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
            p.addOptional('xk_min', xk_min_default, validScalarNum);
            p.addOptional('xk_max', xk_max_default, validScalarNum);
            p.addOptional('xq_min', xq_min_default, validScalarNum);
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

            % non-dimensionalize inputs
            this.rho0        = p.Results.rho0 / this.rho_ref;
            this.p0          = p.Results.p0 ...
                             / (this.rho_ref * this.a_ref^2);
            this.a0        = sqrt( this.gamma * this.p0 / this.rho0 );
            this.q_min     = p.Results.q_min / this.a_ref;

            % already non-dimensional
            this.k_min     = p.Results.k_min;
            this.k_max     = p.Results.k_max;

            if ~isnan(p.Results.xk_min)
                this.xk_min = p.Results.xk_min;
                [this.k_min,~] = this.q_theta_from_x_y(p.Results.xk_min,0);
            else
                [this.xk_min,~] = this.xy_from_phi_psi(0.0,1./this.k_min);
            end

            if ~isnan(p.Results.xk_max)
                this.xk_max = p.Results.xk_max;
                [this.k_max,~] = this.q_theta_from_x_y(this.xk_max,0);
            else
                [this.xk_max,~] = this.xy_from_phi_psi(0.0,1./this.k_max);
            end

            if this.k_min >= this.k_max
                error('k_min must be less than k_max');
            end
            if this.xk_min >= this.xk_max
                error('xk_min must be less than xk_max');
            end

            psi = 1/this.k_max;
            if ~isnan(p.Results.xq_min)
                this.xq_min = p.Results.xq_min;
                if this.xq_min <= this.xk_max
                    error('xq_min must be greater than xk_max')
                end
                [this.q_min,~] = this.q_theta_from_x_psi(this.xq_min,psi);
            else
                if this.q_min >= this.k_min
                    error('q_min must be less than k_min')
                end
                [this.xq_min,~] = this.xy_from_q_psi(this.q_min,psi);
            end
            
        end
        function grid = make_grid_phi_psi(this,N_i,N_j)
            grid = struct();
            grid.imax = N_i;
            grid.jmax = N_j;
            % psi = linspace(1./this.k_max,1./this.k_min,N_j);
            psi = 1./linspace(this.k_max,this.k_min,N_j);
            theta0 = asin(psi(1)*this.q_min);
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
                % theta_max = pi - asin(psi_0(j)*this.q_min);
                theta_max = pi/2;
                theta = linspace(theta_max,theta_min,N_i);
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
        % function s = streamline_arc_length_check_1(this,psi,q0,q1,Nq)
        %     if nargin < 5
        %         Nq = 101;
        %     end
        %     q = linspace(q0,q1,Nq);
        %     [x,y] = this.xy_from_q_psi(q,psi);
        %     s = sum(sqrt(diff(x).^2+diff(y).^2));
        % end
        % function s = streamline_arc_length_check_2(this,psi,q0,q1)
        %     s = integral(@(q)diff_streamline_arc_length(this,psi,q),...
        %                  q0,q1,"AbsTol",1e-15,'RelTol',1e-14);
        %     function ds = diff_streamline_arc_length(this,psi,q)
        %         h = 1e-20;
        %         [dx,dy] = this.xy_from_q_psi(q+1i*h,psi);
        %         dx = imag(dx)/h;
        %         dy = imag(dy)/h;
        %         ds = sqrt(dx.^2+dy.^2);
        %     end
        % end
        % function s = potential_line_arc_length_check_1(this,phi,k0,k1,Nk)
        %     if nargin < 5
        %         Nk = 101;
        %     end
        %     psi = 1./linspace(1/k0,1/k1,Nk);
        %     [x,y] = this.xy_from_phi_psi(phi,psi);
        %     s = sum(sqrt(diff(x).^2+diff(y).^2));
        % end
        % function s = potential_line_arc_length_check_2(this,phi,k0,k1)
        %     s = integral(@(k)diff_potential_line_arc_length(this,phi,k),...
        %         k0,k1,"AbsTol",1e-15,'RelTol',1e-14);
        %     function ds = diff_potential_line_arc_length(this,phi,k)
        %         h = 1e-20;
        %         [dx,dy] = this.xy_from_phi_psi(phi,1./k+1i*h);
        %         dx = imag(dx)/h;
        %         dy = imag(dy)/h;
        %         ds = sqrt(dx.^2+dy.^2);
        %     end
        % end
        
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
        function [q,theta] = q_theta_from_x_psi(this,x,psi)
            opts =  optimset('TolFun',1e-15,'TolX',1e-17);
            q = arrayfun( @(x,k) fzero( @(q)obj_fun(this,x,psi,q),...
                [1e-6,this.k_min],opts),x,psi);
            psi = this.psi_from_x_q(x,q);
            arg = psi.*q;
            arg(arg<-1) = -1;
            arg(arg>1)  =  1;
            theta = asin(arg);
            function f = obj_fun(this,x,psi,q)
                J = this.J_value(q);
                r = this.rho_value(q);
                y = psi./(r.*q) .* sqrt( max(1.0 - (q.*psi).^2,0.0));
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
            M = q.^2./(1-0.5*(this.gamma-1)*q.^2);
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
        function phi = phi_from_q_psi(this,q,psi)
            theta = asin(psi.*q);
            rho = this.rho_value(q);
            phi = cos(theta)./(q.*rho);
        end
    end
end