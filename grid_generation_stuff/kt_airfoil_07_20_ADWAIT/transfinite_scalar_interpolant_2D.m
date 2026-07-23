classdef transfinite_scalar_interpolant_2D
    properties
        G(4,1)  cell
    end
    methods
%% constructor
        function this = transfinite_scalar_interpolant_2D(varargin)
            p = inputParser;
            validReal        = @(x) isreal(x) && ~isnan(x) && ~isinf(x);
            validRealPoints  = @(x) all(arrayfun(@(x)validReal(x),x),"all");
            valid_fcn_handle = @(x) isa(x,'function_handle');
            valid_boundary_input = @(x) valid_fcn_handle(x) || validRealPoints(x);
            addOptional(p,'G1' ,@(xi)                xi, valid_boundary_input);
            addOptional(p,'G2' ,@(eta)              eta, valid_boundary_input);
            addOptional(p,'G3' ,@(xi)                xi, valid_boundary_input);
            addOptional(p,'G4' ,@(eta)              eta, valid_boundary_input);
            parse(p,varargin{:});
            results = structfun(@(x)transfinite_scalar_interpolant_2D.setup_boundaries(x),p.Results,'UniformOutput',false);
            this.G{1,1} = results.G1;
            this.G{2,1} = results.G2;
            this.G{3,1} = results.G3;
            this.G{4,1} = results.G4;
        end
        
        function plot_val(this,n,varargin)
            t = linspace(0,1,n);
            val = this.eval_on_grid(t,t);
            surf(t,t,val)
        end
        function val = eval_on_grid(this,xi_vec,eta_vec)
            val = this.eval(xi_vec(:),eta_vec(:).');
        end
        function val = eval(this,xi,eta)
            val = transfinite_scalar_interpolant_2D.transfinite_eval( xi, eta,     ...
                                                               this.G{1,1}, ...
                                                               this.G{2,1}, ...
                                                               this.G{3,1}, ...
                                                               this.G{4,1} );
        end
    end
    methods (Static)
        function F = setup_boundaries(G)
            if isa(G,'function_handle')
                F = @(x) G( x );
                return
            elseif isscalar(G)
                F = @(x) G*ones(size(x));
            else
                t = linspace(0,1,numel(G));
                pp = csape(t,G);
                F = @(x) fnval(pp,x);
            end
        end
        function val = transfinite_eval(xi,eta,G1,G2,G3,G4)
            x1 = G1(-1);
            x2 = G1( 1);
            x3 = G3( 1);
            x4 = G3(-1);
            val = (1 - xi ).*G4(eta) +  xi.*G2(eta)  ...
              +   (1 - eta).*G1(xi)  + eta.*G3(xi)  ...
              - ( (1 - xi ).* ( (1 - eta).*x1 + eta.*x4 ) ...
              +        xi  .* ( (1 - eta).*x2 + eta.*x3 ) );

        end
    end
end