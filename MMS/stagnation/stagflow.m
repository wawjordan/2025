classdef stagflow
    properties
        rho_ref(1,1)     double   = 1.0
        a_ref(1,1)       double   = 374.165738677
        p_init(1,1)      double   = 100000.0
        gamma(1,1)       double   = 1.4
        
        l_ref(1,1)       double  = 1.0
        l_ref_grid(1,1)  double  = 1.0

        rhoinf(1,1)      double  = 1.0
        mach_max(1,1)    double  = 0.5
        p_min_ratio(1,1) double  = 0.8
        x0(1,1)          double  = 0.0
        y0(1,1)          double  = 0.0
        lx(1,1)          double  = 1.0
        ly(1,1)          double  = 0.5

        pstag(1,1)       double  = 1.0
        vmax(1,1)        double  = 1.0
        vscale(1,1)      double  = 1.0
    end
    methods
%% Constructor
    function this = stagflow(varargin)
        p = inputParser;
        validScalarNum     = @(x) isnumeric(x) && isscalar(x);
        validScalarPosNum  = @(x) validScalarNum(x) && (x > 0);
        addOptional(p,  'rhoinf',      1.0, validScalarPosNum  );
        addOptional(p,  'mach_max',    0.5, validScalarPosNum  )
        addOptional(p,  'p_min_ratio', 0.8, validScalarPosNum  );
        addOptional(p,  'x0',          0.0, validScalarNum     );
        addOptional(p,  'y0',          0.0, validScalarNum  );
        addOptional(p,  'lx',          1.0, validScalarNum  );
        addOptional(p,  'ly',          0.5, validScalarPosNum  );
        addOptional(p,  'p_init',      100000.0, validScalarPosNum  );
        addOptional(p,  'gamma',       1.4, validScalarPosNum  );
        addOptional(p,  'rho_ref',     1.0, validScalarPosNum  );
        addOptional(p,  'a_ref',       374.165738677, validScalarPosNum  );
        addOptional(p,  'l_ref',       1.0, validScalarPosNum  );
        addOptional(p,  'l_ref_grid',  1.0, validScalarPosNum  );
        parse(p,varargin{:});
        this.rho_ref    = p.Results.rho_ref;
        this.a_ref      = p.Results.a_ref;
        this.p_init     = p.Results.p_init / (this.rho_ref * this.a_ref^2);
        this.gamma      = p.Results.gamma;
        this.l_ref      = p.Results.l_ref;
        this.l_ref_grid = p.Results.l_ref_grid;
        L = this.l_ref / this.l_ref_grid;

        this.rhoinf      = p.Results.rhoinf / this.rho_ref;
        this.mach_max    = p.Results.mach_max;
        this.p_min_ratio = p.Results.p_min_ratio;
        this.x0          = p.Results.x0 * L;
        this.y0          = p.Results.y0 * L;
        this.lx          = p.Results.lx * L;
        this.ly          = p.Results.ly * L;

        % convenience vars
        g = this.gamma;
        pm = this.p_min_ratio;
        m2 = this.mach_max^2;
        r  = this.rhoinf;
        xm = this.lx;
        ym = this.ly;
        p0 = this.p_init * (1 + 0.5*(g - 1)*m2)^( g/(g-1) );

        % make sure that maximum mach number and minimum pressure
        % constraints are satisfied
        v1 = sqrt( ( m2*g*p0/r ) / ( 1 + 0.5*g*m2 ) );
        v2 = sqrt( 2*p0*(1-pm) / r );                               
        this.vmax   = min(v1,v2);
        this.vscale = this.vmax / sqrt( xm^2 + ym^2 );
        this.pstag  = p0;
    end
    function u = x_velocity(this,x,~)
        u = this.vscale*(this.x0-x);
    end
    function v = y_velocity(this,~,y)
        v = this.vscale*(y-this.y0);
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
        parse(p,this,nx,ny,varargin{:});
        fx = p.Results.fx;
        fy = p.Results.fy;
        grid = struct();
        grid.imax = nx;
        grid.jmax = ny;
        % x = linspace(this.x0,this.x0+this.lx,nx);
        % y = linspace(this.y0-this.ly,this.y0+this.ly,ny);
        x_ = [this.x0,this.x0+this.lx];
        y_ = [this.y0-this.ly,this.y0+this.ly];
        % x = linspace(min(x_),max(x_),nx);
        % y = linspace(min(y_),max(y_),ny);
        x = fx(min(x_),max(x_),nx);
        y = fy(min(y_),max(y_),ny);
        [grid.x,grid.y] = ndgrid(x,y);
    end
    end
    methods (Static)
        % f = @(x,s,k) (1/(1-(2*s/k)*tanh(k/2)))*( x - (s/k)*(tanh(k*(x-0.5)) + tanh(k/2) ) );
    end
end