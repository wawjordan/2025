classdef ogive_geometry_old
    properties
        L0(1,1)      double = 0.5
        Rn(1,1)      double = 12.5/1000
        Rb(1,1)      double = 45.0/1000
        rho_ogv(1,1) double = 0.0 %2.800277777777778
        x0(1,1)      double = 0.0 %0.075559839317721
        xa(1,1)      double = 0.0 %0.063059839317721
        xt(1,1)      double = 0.0 %0.073656710371330
        yt(1,1)      double = 0.0 %0.012354274572528
        L(1,1)       double = 0.0 %0.436940160682279
        t_circ0(1,1) double = 0.0
        t_circ1(1,1) double = 0.0 %1.417951591790907
        t_ogv0(1,1)  double = 0.0 %
        t_ogv1(1,1)  double = 0.0 %
    end
    methods
        function this = ogive_geometry_old(varargin)
            p = inputParser;
            validScalarNum    = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
            addOptional(p,'L0',     0.50000,validScalarNum);
            addOptional(p,'Rn'     ,0.01250,validScalarPosNum);
            addOptional(p,'Rb'     ,0.04500,validScalarPosNum);
            parse(p,varargin{:});
            this.L0      = p.Results.L0;
            
            this.Rn      = p.Results.Rn;
            this.Rb      = p.Results.Rb;
            this.rho_ogv = (this.Rb^2+this.L0^2)/(2*this.Rb);
            % this.x0      = this.L0 - sqrt( (this.rho_ogv-this.Rn)^2 ...
            %                              - (this.Rb-this.rho_ogv)^2 );
            this.x0      = this.L0 - sqrt( (this.Rb-this.Rn) ...
                                          *(this.L0^2/this.Rb - this.Rn) );
            this.xt      = this.x0 - this.Rn ...
                                   * sqrt( 1 - (( this.L0^2 - this.Rb^2 ) ...
                                               /( this.L0^2 + this.Rb^2  ...
                                              - 2*this.Rb*this.Rn ))^2 );
            this.yt      = ( this.Rn * ( this.L0^2 - this.Rb^2 ) ) ...
                         ./( this.L0^2 + this.Rb^2 - 2*this.Rb*this.Rn );
            this.xa      = this.x0 - this.Rn;
            this.L       = this.L0 - this.xa;
            this.t_circ0  = 0;
            this.t_circ1  = pi - atan2( this.yt, this.xt - this.x0 );
            xc = this.L0;
            yc = (this.Rb^2-this.L0^2)/(2*this.Rb);
            this.t_ogv0   = pi - atan2( this.yt - yc, this.xt - xc );
            this.t_ogv1   = pi/2;
        end
        function y = y_from_x(this,x)
            y = ogive_geometry_old.y_val( this.rho_ogv, ...
                                      this.L0,      ...
                                      this.Rb,      ...
                                      this.Rn,      ...
                                      this.x0,      ...
                                      this.xt,      ...
                                      x );
        end
        function dydx = diff_y_from_x( this, x )
            dydx = ogive_geometry_old.dydx_val( this.rho_ogv, ...
                                            this.L0,      ...
                                            this.Rn,      ...
                                            this.x0,      ...
                                            this.xt,      ...
                                            x );
        end
        function h = plot_ogive( this, varargin )
            h = fplot(@(x) this.y_from_x(x),[this.xa,this.L0],varargin{:});
            xlim([0,this.L0])
            ylim
            axis equal
        end
        function xs = arc_length_param(this,N,F)
            if ( nargin<3)
                F = @(t) t;
            end
            options = optimset('TolFun',1e-15,'TolX',1e-17);
            xmap = @(t) this.xa + this.L * t;
            L_total = this.arc_length( xmap(0), xmap(1) );
            dL = 1 / L_total;
            t  = linspace(0,1,N);
            ts = zeros(N,1);
            ts(1) = 1e-12;
            ts(N) = 1;
            for i = 2:N-1
                ftmp = @(tsi) F( t(i) ) - F( t(i-1) ) ...
                            - dL * this.arc_length( xmap(ts(i-1)), xmap(tsi) );
                ts(i) = fzero( @(tsi)ftmp(tsi),[ts(i-1),1],options);
            end
            ts(1) = 0;
            xs = xmap(ts);
        end
        function S = arc_length(this,xa,xb)
            S = integral( @(t) diff_arc_length(this,t), xa, xb, ...
                          "AbsTol", 1e-15, 'RelTol', 1e-14 );
            function dS = diff_arc_length(this,x)
                dS = sqrt(1 + this.diff_y_from_x(x).^2);
            end
        end
    end
    methods (Static)
        % function xt = get_xt(rho_ogv,L0,Rb,Rn,x0)
        %     options = optimset('TolFun',1e-15,'TolX',1e-17);
        %     xt = fminbnd( @(x) ( ogive_geometry.y_circ( Rn, x0, x ) ...
        %                        - ogive_geometry.y_ogive( rho_ogv, L0, Rb, x )...
        %                        ).^2, x0-Rn, x0+Rn, options );
        % end
        function y = y_val(rho_ogv,L0,Rb,Rn,x0,xt,x)
            y = zeros(size(x));
            mask = (x<=xt);
            y(mask)  = ogive_geometry_old.y_circ(Rn,x0,x(mask));
            y(~mask) = ogive_geometry_old.y_ogive(rho_ogv,L0,Rb,x(~mask));
        end
        function dydx = dydx_val(rho_ogv,L0,Rn,x0,xt,x)
            dydx = zeros(size(x));
            mask = (x<=xt);
            dydx(mask)  = ogive_geometry_old.dydx_circ(Rn,x0,x(mask));
            dydx(~mask) = ogive_geometry_old.dydx_ogive(rho_ogv,L0,x(~mask));
        end
        function y = y_ogive(rho_ogv,L0,Rb,x)
            y = sqrt(rho_ogv^2 - (L0-x).^2) + Rb - rho_ogv;
        end
        function y = y_circ(Rn,x0,x)
            y = sqrt( Rn^2 - (x-x0).^2 );
        end
        function dydx = dydx_ogive(rho_ogv,L0,x)
            dydx = (L0-x)./sqrt(rho_ogv^2 - (L0 - x).^2);
        end
        function dydx = dydx_circ(Rn,x0,x)
            dydx = -(x-x0)./sqrt(Rn^2- (x-x0).^2);
        end
    end
end