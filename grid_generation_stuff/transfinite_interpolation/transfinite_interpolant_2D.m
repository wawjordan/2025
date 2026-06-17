classdef transfinite_interpolant_2D
    properties
        G(4,2)  cell
    end
    methods
%% constructor
        function this = transfinite_interpolant_2D(varargin)
            p = inputParser;
            % validScalarNum = @(x) isnumeric(x) && isscalar(x);
            % validScalarPosNum = @(x) validScalarNum(x) && (x > 0);
            % validScalarPosInt = @(x) mod(x,1)<10*eps(1) && isscalar(x) && (x > 0);
            % validReal        = @(x) isreal(x) && ~isnan(x) && ~isinf(x);
            % validRealPoints  = @(x) all(arrayfun(@(x)validReal(x),x),"all") && isvector(x);
            % validVec         = @(x) validRealPoints(x) && numel(x)==2;
            valid_fcn_handle = @(x) isa(x,'function_handle');
            % addRequired(p,'airfoil_x',validRealPoints);
            % addRequired(p,'airfoil_y',validRealPoints);
            % addRequired(p,'boundary_distance',validScalarPosNum);
            % addRequired(p,'num_wake_pts',validScalarPosInt);
            % addOptional(p,'TE_loc',[1;0],validVec);
            % addOptional(p,'TE_slope',0,validReal);
            % eta min boundary
            addOptional(p,'G1x' ,@(xi)                xi, valid_fcn_handle);
            addOptional(p,'G1y' ,@(xi)  -ones(size(xi)) , valid_fcn_handle);
            addOptional(p,'G2x' ,@(eta)  ones(size(eta)), valid_fcn_handle);
            addOptional(p,'G2y' ,@(eta)              eta, valid_fcn_handle);
            addOptional(p,'G3x' ,@(xi)                xi, valid_fcn_handle);
            addOptional(p,'G3y' ,@(xi)   ones(size(xi)) , valid_fcn_handle);
            addOptional(p,'G4x' ,@(eta) -ones(size(eta)), valid_fcn_handle);
            addOptional(p,'G4y' ,@(eta)              eta, valid_fcn_handle);
            parse(p,varargin{:});
            this.G{1,1} = p.Results.G1x;
            this.G{1,2} = p.Results.G1y;
            this.G{2,1} = p.Results.G2x;
            this.G{2,2} = p.Results.G2y;
            this.G{3,1} = p.Results.G3x;
            this.G{3,2} = p.Results.G3y;
            this.G{4,1} = p.Results.G4x;
            this.G{4,2} = p.Results.G4y;
        end
        function plot_boundaries(this,n,varargin)
            hold on
            t = linspace(-1,1,n);
            x = this.eval_x(t,-1);
            y = this.eval_y(t,-1);
            plot(x,y,varargin{:})

            x = this.eval_x(1,t);
            y = this.eval_y(1,t);
            plot(x,y,varargin{:})

            x = this.eval_x(t,1);
            y = this.eval_y(t,1);
            plot(x,y,varargin{:})

            x = this.eval_x(-1,t);
            y = this.eval_y(-1,t);
            plot(x,y,varargin{:})
        end
        function val = eval_x(this,xi,eta)
            val = transfinite_interpolant_2D.transfinite_eval( xi, eta,     ...
                                                               this.G{1,1}, ...
                                                               this.G{2,1}, ...
                                                               this.G{3,1}, ...
                                                               this.G{4,1} );
        end
        function val = eval_y(this,xi,eta)
            val = transfinite_interpolant_2D.transfinite_eval( xi, eta,     ...
                                                               this.G{1,2}, ...
                                                               this.G{2,2}, ...
                                                               this.G{3,2}, ...
                                                               this.G{4,2} );
        end

    end
    methods (Static)
        function val = transfinite_eval(xi,eta,G1,G2,G3,G4)
            x1 = G1(-1);
            x2 = G1( 1);
            x3 = G3( 1);
            x4 = G3(-1);
            val = 0.5 * ( (1 - xi ).*G4(eta) + (1 +  xi).*G2(eta)   ...
                      +   (1 - eta).*G1(xi)  + (1 + eta).*G3(xi)  ) ...
               - 0.25 * ( (1 - xi ).* ( (1 - eta).*x1 + (1 + eta).*x4 ) ...
                      +   (1 + xi ).* ( (1 - eta).*x2 + (1 + eta).*x3 ) );

        end
    end
end