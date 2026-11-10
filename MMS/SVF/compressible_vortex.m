classdef compressible_vortex
    properties
        theta_0(1,1)     double = 90.0
        theta_1(1,1)     double = 0.0
        r_inner(1,1)     double = 2.0
        r_outer(1,1)     double = 3.0
        eta(1,1)         double = 1.5
        rho_inner(1,1)   double = 1.0
        p_inner(1,1)     double = 0.8611583247416177E+05
        M_inner(1,1)     double = 2.0
        v_inner(1,1)     double = 1.0
        gamma(1,1)       double = 1.4
        gm1x2(1,1)       double = 0.2
        rho_ref(1,1)     double = 1.0
        a_ref(1,1)       double = 340.0
        l_ref(1,1)       double = 1.0
        l_ref_grid(1,1)  double = 1.0
    end
    methods
%% Constructor
        function this = compressible_vortex(varargin)
            validScalarNum    = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
            p = inputParser;
            p.addOptional('r_inner',2.0,validScalarPosNum);
            p.addOptional('r_outer',3.0,validScalarPosNum);
            p.addOptional('theta_0',90.0,validScalarNum);
            p.addOptional('theta_1', 0.0,validScalarNum);
            p.addOptional('M_inner',  2.0,validScalarPosNum)
            p.addOptional('rho_inner',1.0,validScalarPosNum)
            p.addOptional('p_inner',-999,validScalarPosNum)
            p.addOptional('gamma',1.4,validScalarPosNum)
            p.addOptional('rho_ref',1.0,validScalarPosNum)
            p.addOptional('a_ref',340.0,validScalarPosNum)
            p.addOptional('l_ref',1.0,validScalarPosNum)
            p.addOptional('l_ref_grid',1.0,validScalarPosNum)
            parse(p,varargin{:});
            this.l_ref       = p.Results.l_ref;
            this.l_ref_grid  = p.Results.l_ref_grid;
            L = this.l_ref / this.l_ref_grid;
            this.theta_0     = deg2rad( p.Results.theta_0 );
            this.theta_1     = deg2rad( p.Results.theta_1 );
            this.r_inner     = p.Results.r_inner / L;
            this.r_outer     = p.Results.r_outer / L;
            this.eta         = this.r_outer / this.r_inner;
            this.rho_ref     = p.Results.rho_ref;
            this.a_ref       = p.Results.a_ref;
            this.gamma       = p.Results.gamma;
            this.gm1x2       = (this.gamma - 1)/2;
            this.rho_inner   = p.Results.rho_inner / this.rho_ref;
            if p.Results.p_inner < 0
                this.p_inner = 1/this.gamma;
            else
                this.p_inner = p.Results.p_inner ...
                             / (this.rho_ref * this.a_ref^2);
            end
            this.M_inner = p.Results.M_inner;
            this.v_inner =  this.r_inner * this.M_inner ...
                                         * this.rho_inner^this.gm1x2;
        end
        function grid = make_grid(this,N_theta,N_r)
            grid = struct();
            grid.imax = N_theta;
            grid.jmax = N_r;
            theta = linspace( this.theta_0, this.theta_1, N_theta );
            r     = compressible_vortex.r_distribution( this.r_inner, ...
                                                        this.r_outer, ...
                                                        linspace(0,1,N_r));
            [THETA,R] = ndgrid(theta,r);
            grid.x = R.*cos(THETA);
            grid.y = R.*sin(THETA);
        end
        function rho = density(this,x,y)
            r2 = x.^2 + y.^2;
            rho = this.r_inner * ( 1 + this.gm1x2*this.M_inner^2 ...
                                      *( 1 - (this.r_inner^2./r2) ) ...
                                      ).^(1/(this.gamma-1));
        end
        function u = x_velocity(this,x,y)
            r2 = x.^2 + y.^2;
            u  = this.v_inner*y./r2;
        end
        function v = y_velocity(this,x,y)
            r2 = x.^2 + y.^2;
            v  =-this.v_inner*x./r2;
        end
        function p = pressure(this,x,y)
            p = this.density(x,y).^(this.gamma)/this.gamma;
        end
    end
    methods (Static)
        function AR = calc_AR_simple(r_0,r_1,theta_0,theta_1,N_r,N_theta)
            eta = r_1/r_0;
            deta = eta^(1/(N_r-1));
            dtheta = abs(theta_1 - theta_0)/(N_theta-1);
            AR   = tan(dtheta/2) .* (deta+1)./(deta-1);
        end
        function r = r_distribution(r_0,r_1,tau)
            eta = r_1 / r_0;
            r = (r_1-r_0)*( eta.^tau - 1 )/( eta - 1 )  + r_0;
        end
        % function [AR,deta] = calc_AR( dtheta_0, dtheta_mid )
        %     deta = (1+tan(dtheta_0/2))./(1-tan(dtheta_0/2));
        %     AR   = tan(dtheta_mid/2) .* (deta+1)./(deta-1);
        % end
    end
end