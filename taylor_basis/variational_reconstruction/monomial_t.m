classdef monomial_t

    properties
        exponents  (:,1) {mustBeInteger}        % monomial term exponents
    end

    methods

        function this = monomial_t(exponents)
            if nargin > 0
                this.exponents = exponents(:);
            end
        end

        function [val,coef] = eval(this,x)
            [val,coef] = monomial_t.static_eval(numel(this.exponents),this.exponents,x);
        end

        function [dval,coef] = deval(this,x,order)
            [dval,coef] = monomial_t.static_eval_derivative( ...
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

        function [dval,coef] = static_eval_derivative(n_dim,exponents,x,order)
            diff_order = exponents(:) - order(:);
            if any( diff_order < 0 )
                coef = 1;
                dval = 0.0;
                return;
            end
            [dval,coef] = monomial_t.static_eval(n_dim,diff_order,x);
        end
    end

end