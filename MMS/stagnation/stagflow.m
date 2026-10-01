classdef stagflow
    properties
        l(1,1)         double  = 1.0
        x0(1,1)        double  = 0.0
        y0(1,1)        double  = 0.0
        lx(1,1)        double  = 1.0
        ly(1,1)        double  = 0.5
        gamma(1,1)     double  = 1.4
        pref(1,1)      double  = 1.0
        rhoref(1,1)    double  = 1.0
        aref(1,1)      double  = 1.0
        vref(1,1)      double  = 1.0
        rhoinf(1,1)    double  = 1.0
        mach_max(1,1)  double  = 0.5
        pstag(1,1)     double  = 1.0
        pmin(1,1)      double  = 0.5
        vmax(1,1)      double  = 1.0
    end  
    methods
%% Constructor
    function this = stagflow(varargin)
            p = inputParser;
            validScalarNum     = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum  = @(x) validScalarNum(x) && (x > 0);
            addOptional(p,  'l',          1.0, validScalarPosNum  );
            addOptional(p,  'x0',         0.0, validScalarNum     );
            addOptional(p,  'y0',         0.0, validScalarNum  );
            addOptional(p,  'lx',         1.0, validScalarNum  );
            addOptional(p,  'ly',         0.5, validScalarPosNum  );
            addOptional(p,  'mach_max',   0.9, validScalarPosNum  );
            addOptional(p,  'aref',       1.0, validScalarPosNum  );
            addOptional(p,  'pstag',      1.0, validScalarPosNum  );
            addOptional(p,  'pmin',       0.9, validScalarPosNum  );
            addOptional(p,  'rhoref',     1.0, validScalarPosNum  );
            addOptional(p,  'rhoinf',     1.0, validScalarPosNum  );
            addOptional(p,  'gamma',      1.4, validScalarPosNum  );
            parse(p,varargin{:});
            if ( p.Results.pmin >= p.Results.pstag )
                error('pmin must be less than pstag')
            end
            this.l        = p.Results.l;
            this.x0       = p.Results.x0 / this.l;
            this.y0       = p.Results.y0 / this.l;
            this.lx     = p.Results.lx / this.l;
            this.ly     = p.Results.ly / this.l;
            this.gamma    = p.Results.gamma;
            this.rhoref   = p.Results.rhoref;
            this.rhoinf   = p.Results.rhoinf / this.rhoref;
            this.aref     = p.Results.aref;
            pref_ = (this.rhoref*this.aref^2);
            this.pstag    = p.Results.pstag / pref_;
            this.pmin     = p.Results.pmin  / pref_;
            this.mach_max = p.Results.mach_max;

            % convenience vars
            g = this.gamma;
            p = this.pstag;
            pm = this.pmin;
            m2 = this.mach_max^2;
            r = this.rhoinf;

            % make sure that maximum mach number and minimum pressure
            % constraints are satisfied
            vmax1 = sqrt( g*(p/r) * m2/(1+0.5*g*m2) );
            vmax2 = sqrt( (p-pm)/(0.5*r) );
            this.vmax = min(vmax1,vmax2);
            this.vref = this.vmax/sqrt(this.lx.^2+this.ly^2);
    end
    function u = x_velocity(this,x,~)
        u = this.vref*(this.x0-x);
    end
    function v = y_velocity(this,~,y)
        v = this.vref*(y-this.y0);
    end
    function p = pressure(this,x,y)
        u = this.x_velocity(x,y);
        v = this.y_velocity(x,y);
        p  = this.pstag - 0.5*this.rhoinf*(u.^2 + v.^2);
    end
    function rho = density(this,x,~)
        rho = this.rhoinf*ones(size(x));
    end
    function m   = mach(this,x,y)
        vmag2 = this.x_velocity(x,y).^2 + this.y_velocity(x,y).^2;
        p    = this.pstag - 0.5*this.rhoinf*vmag2;
        m    = sqrt(vmag2)./sqrt( this.gamma*p/this.rhoinf );
    end
    function grid = make_grid(this,nx,ny)
        grid = struct();
        grid.imax = nx;
        grid.jmax = ny;
        x = linspace(this.x0,this.x0+this.lx,nx);
        y = linspace(this.y0-this.ly,this.y0+this.ly,ny);
        [grid.x,grid.y] = ndgrid(x,y);
    end
    end
    methods (Static)
    end
end