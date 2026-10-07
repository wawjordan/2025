classdef stagflow2
    properties
        length(1,1)      double = 1.0
        l(1,1)           double = 0.0
        R(1,1)           double = 1.0
        alpha(1,1)       double = 0.0
        rho_inf(1,1)     double = 1.0
        v_inf(1,1)       double = 68.0
        p_inf(1,1)       double = 100000.0
        gamma(1,1)       double = 1.4
        rho_ref(1,1)     double = 1.0
        a_ref(1,1)       double = 340.0
        l_ref(1,1)       double = 1.0
        l_ref_grid(1,1)  double = 1.0
        % v_y=v_inf ~ 0.980580675690920*l
        % 1 = Re( 1i*( exp(-1i*alpha) - R^2*exp(1i*alpha)./zeta.^2 )/(1 - l^2./zeta.^2) )
        % 1 = Re( 1i*( -1i - 1i*l^2./zeta.^2 )/(1 - l^2./zeta.^2) )
        % 1 = Re( 1i*(-1i)*( 1 + l^2./zeta.^2 )/(1 - l^2./zeta.^2) )
        % 1 = Re( ( 1 + l^2./zeta.^2 )/(1 - l^2./zeta.^2) )
        % zeta.^2 = (zr +1i*zi)(zr +1i*zi) = zr^2 + 2i*zr*zi - zi^2
        % 1 = Re{ [ ( zeta.^2 * (1 + l^2) )/(l^2zeta.^2) ]/[ ( zeta.^2 * (1 - l^2) )/(l^2zeta.^2) ] }
        % 1 = Re{ [ ( zeta.^2 * (1 + l^2) ) ]/[ ( zeta.^2 * (1 - l^2) ) ]
        % (1 - l^2./zeta.^2) = ( 1 + l^2./zeta.^2 )
        % 0 = l^2./zeta.^2
    end
    methods
%% Constructor
        function this = stagflow2(varargin)
            validScalarNum        = @(x) isnumeric(x) && isscalar(x);
            validScalarNonNegNum  = @(x) validScalarNum(x) && (x >= 0);
            validScalarPosNum     = @(x) validScalarNum(x) && (x > 0);
            p = inputParser;
            p.addOptional('length',1.0,validScalarNonNegNum);
            p.addOptional('alpha',0.0,validScalarNum)
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
            % modified for flow impacting perpendicular flat plate
            % alpha = alpha + 90 degrees
            this.alpha       = deg2rad( p.Results.alpha ) + pi/2;
            L = this.l_ref / this.l_ref_grid;
            this.length      = p.Results.length / L;
            this.R           = this.length/2;
            this.l           = this.R;
        end
        function z = z_from_zeta(this,zeta)
            z = stagflow2.joukowsky_transform(this.l,zeta);
        end
        function dzdzeta = diff_z_from_zeta(this,zeta)
            dzdzeta = stagflow2.joukowsky_transform_derivative(this.l,zeta);
        end
        function zeta = zeta_from_z(this,z)
            % rotate the coordinate system clockwise by 90 degrees
            zeta = stagflow2. ...
                   inverse_joukowsky_transform( this.l, -1i*z );
        end
        function F = complex_potential_from_zeta(this,zeta)
            F = stagflow2.cylinder_potential( this.v_inf,       ...
                                              this.R,           ...
                                              this.alpha,       ...
                                              zeta );
        end
        function w = complex_velocity_from_zeta(this,zeta)
            dzdzeta = this.diff_z_from_zeta(zeta);
            w = stagflow2.cylinder_velocity( this.v_inf,       ...
                                             this.R,           ...
                                             this.alpha,       ...
                                             zeta )./dzdzeta;
        end
        function F = complex_potential( this, z )
            F = this.complex_potential_from_zeta( this.zeta_from_z(z) );
        end
        function w = complex_velocity( this, z )
            w = this.complex_velocity_from_zeta( this.zeta_from_z(z) );
            w = 1i*w; % rotate it back (counterclockwise)
        end
        function rho = density( this, x, ~ )
            rho = this.rho_inf*ones(size(x));
        end
        function u = x_velocity( this, x, y )
            u = real( this.complex_velocity( x + 1i*y ) );
        end
        function v = y_velocity( this, x, y )
            v = -imag( this.complex_velocity( x + 1i*y ) );
        end
        function p = pressure( this, x, y )
            w2 = abs(this.complex_velocity( x + 1i*y )).^2;
            p  = this.p_inf + 0.5*this.rho_inf*( this.v_inf^2 - w2 );
        end
        function h = plot_plate(this,varargin)
            h = plot([0,0],[-this.length,this.length],varargin{:});
        end
        function grid = make_grid(this,nx,ny,varargin)
            p = inputParser;
            validScalarNum     = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum  = @(x) validScalarNum(x) && (x > 0);
            p.addRequired('this');
            p.addRequired('nx',validScalarPosNum);
            p.addRequired('ny',validScalarPosNum);
            fdefault = @(a,b,n) linspace(a,b,n);
            p.addOptional('fx',fdefault,@(x)isa(x,"function_handle"))
            p.addOptional('fy',fdefault,@(x)isa(x,"function_handle"))
            p.addOptional('x0',-1,validScalarNum)
            p.addOptional('x1', 0,validScalarNum)
            p.addOptional('y0',-0.5,validScalarNum)
            p.addOptional('y1',-0.5,validScalarNum)
            parse(p,this,nx,ny,varargin{:});
            fx = p.Results.fx;
            fy = p.Results.fy;
            x0 = p.Results.x0;
            x1 = p.Results.x1;
            y0 = p.Results.y0;
            y1 = p.Results.y1;
            grid = struct();
            grid.imax = nx;
            grid.jmax = ny;
            x_ = [x0,x1];
            y_ = [y0,y1];
            x = fx(min(x_),max(x_),nx);
            y = fy(min(y_),max(y_),ny);
            [grid.x,grid.y] = ndgrid(x,y);
        end
    end
    methods (Static)
        function F = cylinder_potential(v_inf,R,alpha,zeta)
            F = v_inf * ( zeta*exp(-1i*alpha) + R^2*exp(1i*alpha)./zeta );
        end
        function w = cylinder_velocity(v_inf,R,alpha,zeta)
            w = v_inf * ( exp(-1i*alpha) - R^2*exp(1i*alpha)./zeta.^2 );
        end
        function z = joukowsky_transform(l,zeta)
            tol  = 1.0e-12;
            z    = zeros(size(zeta));
            mask = abs(zeta)>tol;
            z(mask)  = zeta(mask) + l^2./zeta(mask);
        end
        function dzdzeta = joukowsky_transform_derivative(l,zeta)
            tol     = 1.0e-12;
            dzdzeta = ones(size(zeta));
            mask    = abs(zeta)>tol;
            dzdzeta(mask) = 1 - l^2./zeta(mask).^2;
        end
        function zeta = inverse_joukowsky_transform(l,z)
            tol  = 1.0e-12;
            if l>tol
                zp   = sqrt(z+2*l);
                zm   = sqrt(z-2*l);
                zeta = l*( (zp+zm)./(zp-zm) );
            else
                zeta = z;
            end
        end
    end
end