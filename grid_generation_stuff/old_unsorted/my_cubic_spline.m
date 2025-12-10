classdef my_cubic_spline
    properties
        x       (:,1) {mustBeNumeric}
        y       (:,1) {mustBeNumeric}
        N       (1,1) {mustBeInteger}
        y1      (1,1) {mustBeNumeric}
        yN      (1,1) {mustBeNumeric}
        d2y     (:,1) {mustBeNumeric}
    end
    methods
        function this = my_cubic_spline(x,y,y1,yN)
            this.N = length(x);
            if (x(1)>x(this.N))
                x = flip(x);
                y = flip(y);
                y1 = yN;
                yN = y1;
            end
            this.x = x;
            this.y = y;
            this.d2y = this.myspline(x,y,y1,yN);
        end
        % function y = eval(this,x)
        %     y = splint(this.x,this.y,this.d2y,x);
        % end
        function y = eval(this,x)
            y = arrayfun(@(x) this.splint_eval(this.x,this.y,this.d2y,x), x );
        end
        function y = deval(this,x)
            y = arrayfun(@(x) this.splint_deval(this.x,this.y,this.d2y,x), x );
        end
    end
    methods(Static)
        function y2 = myspline(x,y,yp1,ypn)

            n = length(x);
            u = zeros(n,1);
            y2 = zeros(n,1);
            LARGE = 0.99e30;

            if ( yp1 > LARGE )
                y2(1) = 0;
                u(1) = 0;
            else
                y2(1) = -0.5;
                u(1) = (3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1);
            end
            for i = 2:n-1
                sig = (x(i)-x(i-1))/(x(i+1)-x(i-1));
                p = sig*y2(i-1)+2;
                y2(i) = (sig-1)/p;
                u(i) = (6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))...
                    /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p;
            end
            if ( ypn > LARGE )
                qn = 0;
                un = 0;
            else
                qn = 0.5;
                un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)));
            end
            y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1);
            for k = n-1:-1:1
                y2(k) = y2(k)*y2(k+1) + u(k);
            end

        end

        % function y = splint(xa,ya,y2a,x)
        % 
        %     n = length(xa);
        % 
        %     klo = 1;
        %     khi = n;
        % 
        %     while (khi-klo > 1)
        %         k = uint16((khi+klo)/2);
        %         if (xa(k) > x)
        %             khi = k;
        %         else
        %             klo = k;
        %         end
        %     end
        % 
        %     h = xa(khi) - xa(klo);
        %     if (h == 0)
        %         error('bad xa input in splint')
        %     end
        %     a = (xa(khi) - x)/h;
        %     b = (x - xa(klo))/h;
        %     y=a*ya(klo)+b*ya(khi)+...
        %         ((a^3-a)*y2a(klo)+(b^3-b)*y2a(khi))*(h^2)/6;
        % end

        function y = splint_eval(xa,ya,y2a,x)
            [a,b,h,khi,klo] = my_cubic_spline.splint_core(xa,x);
            y = a*ya(klo) + b*ya(khi) + ( (a^3-a)*y2a(klo) + ...
                                          (b^3-b)*y2a(khi) )*(h^2)/6;
        end

        function dy = splint_deval(xa,ya,y2a,x)
            [a,b,h,khi,klo] = my_cubic_spline.splint_core(xa,x);
            dy = ( ya(khi) - ya(klo) )/h - ( (3*a^2 - 1)*y2a(klo) ...
                                           - (3*b^2 - 1)*y2a(khi) )*h/6;
        end

        function [a,b,h,khi,klo] = splint_core(xa,x)
            n = length(xa);

            klo = 1;
            khi = n;

            while (khi-klo > 1)
                k = uint16((khi+klo)/2);
                if (xa(k) > x)
                    khi = k;
                else
                    klo = k;
                end
            end

            h = xa(khi) - xa(klo);
            if (h == 0)
                error('bad xa input in splint_core')
            end
            a = (xa(khi) - x)/h;
            b = (x - xa(klo))/h;
        end

    end

end