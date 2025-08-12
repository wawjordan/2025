classdef monomial_t6

    properties
        exponents  (:,1) {mustBeInteger}        % monomial term exponents
    end

    methods

        function this = monomial_t6(exponents)
            if nargin > 0
                this.exponents = exponents(:);
            end
        end

        function [val,coef] = eval(this,x)
            [val,coef] = monomial_t6.static_eval(numel(this.exponents),this.exponents,x);
        end

        function [dval,dcoef,coef] = deval(this,x,order)
            [dval,dcoef,coef] = monomial_t6.static_eval_derivative( ...
                                                       numel(this.exponents),...
                                                             this.exponents, ...
                                                             x, order );
        end

    end
    methods (Static)
        function [val,coef] = static_eval(n_dim,exponents,x)
            val = 1.0;
            coef = 1;
            for d = 1:n_dim
                for i = exponents(d):-1:1
                    val = val*x(d);
                    coef = coef * i;
                end
            end
        end

        function [dval,dcoef,coef] = static_eval_derivative(n_dim,exponents,x,order)
            diff_order = exponents(:) - order(:);
            if any( diff_order < 0 )
                dval  = 0.0;
                dcoef = 1;
                coef  = 1;
                return;
            else
                dval  = 1.0;
                dcoef = 1;
                coef  = 1;
                for d = 1:n_dim
                    for i = exponents(d):-1:diff_order(d)+1
                        dcoef = dcoef * i;
                    end
                    for i = diff_order(d):-1:1
                        coef = coef * i;
                        dval = dval*x(d);
                    end
                end
            end
        end
    end

end