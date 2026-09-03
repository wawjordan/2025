classdef ogive_geometry
    properties
        L0(1,1)      double = 0.5
        Rn(1,1)      double = 12.5/1000
        Rb(1,1)      double = 45.0/1000
        rho_ogv(1,1) double = 0.0 %2.800277777777778
        x0(1,1)      double = 0.0 %0.075559839317721
        xa(1,1)      double = 0.0 %0.063059839317721
        xt(1,1)      double = 0.0 %0.073656710371330
        yt(1,1)      double = 0.0 %0.012354274572528
        xoc(1,1)     double = 0.0 %0.500000000000000
        yoc(1,1)     double = 0.0 %-2.755277777777778
        xcc(1,1)     double = 0.0 %0.075559839317721
        ycc(1,1)     double = 0.0 %0.000000000000000
        L(1,1)       double = 0.0 %0.436940160682279
        t_circ0(1,1) double = 0.0 %0.000000000000000
        t_circ1(1,1) double = 0.0 %1.417951591790907
        t_ogv0(1,1)  double = 0.0 %
        t_ogv1(1,1)  double = 0.0 %
    end
    methods
        function this = ogive_geometry(varargin)
            p = inputParser;
            validScalarNum    = @(x) isnumeric(x) && isscalar(x);
            validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
            addOptional(p,'L0',     0.50000,validScalarNum);
            addOptional(p,'Rb'     ,0.04500,validScalarPosNum);
            addOptional(p,'Rn'     ,0.01250,validScalarPosNum);
            parse(p,varargin{:});
            this.L0      = p.Results.L0;
            this.Rb      = p.Results.Rb;
            this.Rn      = p.Results.Rn;
            this.rho_ogv = (this.Rb^2+this.L0^2)/(2*this.Rb);
            this.x0      = this.L0 - sqrt( (this.Rb-this.Rn) ...
                                          *(this.L0^2/this.Rb - this.Rn) );
            % x0      = this.L0 - sqrt( (this.rho_ogv-this.Rn)^2 ...
            %                              - (this.Rb-this.rho_ogv)^2 );
            this.xt      = this.x0 - this.Rn ...
                                   * sqrt( 1 - (( this.L0^2 - this.Rb^2 ) ...
                                               /( this.L0^2 + this.Rb^2  ...
                                              - 2*this.Rb*this.Rn ))^2 );
            this.yt      = ( this.Rn * ( this.L0^2 - this.Rb^2 ) ) ...
                         ./( this.L0^2 + this.Rb^2 - 2*this.Rb*this.Rn );
            this.xa      = this.x0 - this.Rn;
            this.L       = this.L0 - this.xa;
            this.xcc     = this.x0;
            this.ycc     = 0;
            this.xoc     = this.L0;
            this.yoc     = (this.Rb^2-this.L0^2)/(2*this.Rb);
            this.t_circ0 = 0;
            this.t_circ1 = pi - atan2( this.yt - this.ycc, this.xt - this.xcc );
            this.t_ogv0  = pi - atan2( this.yt - this.yoc, this.xt - this.xoc );
            this.t_ogv1  = pi/2;
        end
        function [x,y] = get_coords(this,N_tip,N_ogive,refine)
            if (nargin<4)
                refine=true;
            end
            theta = this.theta_param(N_tip,N_ogive,refine);
            x = this.x_from_theta(theta);
            y = this.y_from_theta(theta);
        end
        function x = x_from_theta( this, theta )
            x        = zeros(size(theta));
            mask     = abs(theta)<=this.t_circ1;
            x(mask)  = this.xcc - this.Rn*cos(theta(mask));
            x(~mask) = this.xoc - this.rho_ogv*cos(theta(~mask));
        end
        function y = y_from_theta( this, theta )
            y        = zeros(size(theta));
            mask     = abs(theta)<=this.t_circ1;
            y(mask)  = this.ycc + this.Rn*sin(theta(mask));
            y(~mask) = this.yoc + this.rho_ogv*sin(theta(~mask));
        end
        function y = y_from_x( this, x )
            y = ogive_geometry.y_val( this.rho_ogv, ...
                                      this.L0,      ...
                                      this.Rb,      ...
                                      this.Rn,      ...
                                      this.x0,      ...
                                      this.xt,      ...
                                      x );
        end
        function h = plot_ogive( this, varargin )
            h = fplot( @(theta) this.x_from_theta(theta), ...
                       @(theta) this.y_from_theta(theta), ...
                       [0,pi/2],varargin{:});
            axis equal
        end
        function theta = theta_param(this,N_tip,N_ogive,refine)
            if (nargin<4)
                refine=true;
            end

            % get arc length of the blunted tip
            s_tip = this.t_circ1*this.Rn;

            % get arc length spacing on the tip
            ds_tip = s_tip/(N_tip-1);

            % get arc length on the rest of the body
            s_ogive = (this.t_ogv1-this.t_ogv0)*this.rho_ogv;

            % get the relative spacing on the body to match the tip
            d_ogive_0 = ds_tip/s_ogive;

            % get total number of points and allocate theta
            N_total = N_tip + N_ogive - 1;
            theta = zeros(N_total,1);

            % theta on the tip is uniformly spaced
            theta(1:N_tip) = linspace(0,this.t_circ1,N_tip);

            % use a stretching function for the remaining points on the body
            s0 = ogive_geometry.vinokur_one_sided_set_spacing(N_ogive,d_ogive_0,refine);

            % compute the stretching function parameter values
            t0 = linspace(0,1,N_ogive);
            t = ogive_geometry.vinokur_one_sided(t0,s0,refine);
            
            % transform these values into the corresponding thetas
            theta_tmp = this.t_ogv0 + (this.t_ogv1-this.t_ogv0)*t;

            % append these thetas
            theta(N_tip+1:end) = theta_tmp(2:end);
        end
    end
    methods (Static)
        function y = y_val(rho_ogv,L0,Rb,Rn,x0,xt,x)
            y = zeros(size(x));
            mask = (x<=xt);
            y(mask)  = ogive_geometry_old.y_circ(Rn,x0,x(mask));
            y(~mask) = ogive_geometry_old.y_ogive(rho_ogv,L0,Rb,x(~mask));
        end
        function y = y_ogive(rho_ogv,L0,Rb,x)
            y = sqrt(rho_ogv^2 - (L0-x).^2) + Rb - rho_ogv;
        end
        function y = y_circ(Rn,x0,x)
            y = sqrt( Rn^2 - (x-x0).^2 );
        end
        function s0 = vinokur_one_sided_set_spacing(N,d0,refine)
            dt = 1/(N-1);
            s =d0/dt;
            s0 = fzero(@(s)obj_fun(N,s,d0,refine),s);
            function f = obj_fun(N,s0,d0,refine)
                d = ogive_geometry.vinokur_one_sided(1/(N-1),s0,refine);
                f = d-d0;
            end
        end
        function xi = vinokur_one_sided(t,s0,refine)
            s0_ = 1/s0;
            [delta,branch] = ogive_geometry.vinokur_get_delta_one_sided(s0_,refine);
            switch branch
                case 0
                    xi  = t.*(1-0.5*(s0_-1)*(1-t).*(2-t));
                case 1
                    xi = 1 + ( tanh(delta*(t-1))/tanh(delta) );
                case 2
                    xi = 1 + ( tan(delta*(t-1))/tan(delta) );
            end
        end
        function [delta,branch] = vinokur_get_delta_one_sided(s0,refine)
            tol = 0.001;
            if abs(s0-1) < tol
                delta  = 0;
                branch = 0;
            elseif s0 > 1
                delta = 0.5*vinokur_sinh_function( s0 );
                if (refine)
                    delta = 0.5*fzero(@(z)sinh(z)./z - s0,2*delta);
                end
                branch = 1;
            elseif s0 < 1
                delta = 0.5*vinokur_sine_function( s0 );
                if (refine)
                    delta = 0.5*fzero(@(z)sin(z)./z - s0,2*delta);
                end
                branch = 2;
            end
            function [ x ] = vinokur_sine_function( y )
                % Approximately solves the function y = sin(x)/x for x
                % if y < 0.26938972
                %     x = sine_function_1(y);
                % else
                %     x = sine_function_2(y);
                % end
                mask = (y < 0.26938972);
                x = zeros(size(y));
                x(mask)  = sine_function_1(y(mask));
                x(~mask) = sine_function_2(y(~mask));
                function x = sine_function_1(y)
                    a =  1;
                    b = -1;
                    c =  1;
                    d = -(1 + pi^2/6);
                    e =  6.794732;
                    f = -13.205501;
                    g =  11.726095;
                    x = pi*(a + y.*(b + y.*(c + y.*(d + y.*(e + y.*(g*y + f))))));
                end

                function x = sine_function_2(y)
                    y1 = 1 - y;
                    a =  1;
                    b =  0.15;
                    c =  0.057321429;
                    d =  0.048774238;
                    e = -0.053337753;
                    f =  0.075845134;
                    x = sqrt(6*y1).*(a + y1.*(b + y1.*(c + y1.*(d + y1.*(f*y1 + e)))));
                end
            end
            function [ x ] = vinokur_sinh_function( y )
                % Approximately solves the function y = sinh(x)/x for x
                % if y < 2.7829681
                %     x = sinh_function_1(y);
                % else
                %     x = sinh_function_2(y);
                % end
                mask = (y < 2.7829681);
                x = zeros(size(y));
                x(mask)  = sinh_function_1(y(mask));
                x(~mask) = sinh_function_2(y(~mask));
                function x = sinh_function_1(y)
                    y1 = y - 1;
                    a =  1;
                    b = -0.15;
                    c =  0.057321429;
                    d = -0.024907295;
                    e =  0.0077424461;
                    f = -0.0010794123;
                    x  = sqrt(6*y1).*(a + y1*(b + y1*(c + y1*(d + y1*(f*y1 + e)))));
                end

                function x = sinh_function_2(y)
                    v = log(y);
                    w = 1./y - 0.028527431;
                    a = -0.02041793;
                    b =  0.24902722;
                    c =  1.9496443;
                    d = -2.6294547;
                    e =  8.56795911;
                    x = v + ( 1 + 1./v ).*log(2*v) + a + w.*(b + w.*(c + w.*(e*w + d)));
                end
            end
        end
    end
end