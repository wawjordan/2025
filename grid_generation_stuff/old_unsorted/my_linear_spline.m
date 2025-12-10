classdef my_linear_spline
    properties
        x       (:,1) {mustBeNumeric}
        y       (:,1) {mustBeNumeric}
        N       (1,1) {mustBeInteger}
        ds      (:,1) {mustBeNumeric}
    end
    methods
        function this = my_linear_spline(x,y)
            N = length(x);
            this.N = N;
            if (x(1)>x(N))
                x = flip(x);
                y = flip(y);
            end
            this.x = x;
            this.y = y;
            this.ds = abs( y(1:N-1) - y(2:N) );
        end
        function y = eval(this,x)
            y = arrayfun(@(x) this.splint_eval(this.x,this.y,x), x );
        end
        function y = deval(this,x)
            y = arrayfun(@(x) this.splint_deval(this.x,this.y,x), x );
        end
        function s = L_approx(this,x1,x2)
            s = this.length_approx(this.x,this.y,this.ds,x1,x2);
        end
    end
    methods(Static)
        function y = splint_eval(xa,ya,x)
            [a,b,~,khi,klo] = my_linear_spline.splint_core(xa,x);
            y = a*ya(klo) + b*ya(khi);
        end

        function dy = splint_deval(xa,ya,x)
            [~,~,h,khi,klo] = my_linear_spline.splint_core(xa,x);
            dy = ( ya(khi) - ya(klo) )/h;
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

        function S = length_approx(xa,ya,ds,x1,x2)
            ds = ds(:);
            flp = false;
            if (x1 > x2)
                xtmp = x2;
                x2 = x1;
                x1 = xtmp;
                flp = true;
            end
            [a,b,~,khi1,klo1] = my_linear_spline.splint_core(xa,x1);
            y1 = a*ya(klo1) + b*ya(khi1);
            ds1 = abs(y1-ya(klo1));
            

            [a,b,~,khi2,klo2] = my_linear_spline.splint_core(xa,x2);
            y2 = a*ya(klo2) + b*ya(khi2);

            ds2 = abs(y2-ya(klo2));
            S = [ds1;ds(khi1:klo2);ds2];
            if flp
                S = flip(S);
            end
        end

    end

end