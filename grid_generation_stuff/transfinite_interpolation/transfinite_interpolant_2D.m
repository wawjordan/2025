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
            validReal        = @(x) isreal(x) && ~isnan(x) && ~isinf(x);
            validRealPoints  = @(x) all(arrayfun(@(x)validReal(x),x),"all");
            validVec         = @(x) validRealPoints(x) && numel(x)==2;
            valid_fcn_handle = @(x) isa(x,'function_handle');
            valid_boundary_input = @(x) valid_fcn_handle(x) || validRealPoints(x);
            addOptional(p,'Gx1' ,@(xi)                xi, valid_boundary_input);
            addOptional(p,'Gy1' ,@(xi)  -ones(size(xi)) , valid_boundary_input);
            addOptional(p,'Gx2' ,@(eta)  ones(size(eta)), valid_boundary_input);
            addOptional(p,'Gy2' ,@(eta)              eta, valid_boundary_input);
            addOptional(p,'Gx3' ,@(xi)                xi, valid_boundary_input);
            addOptional(p,'Gy3' ,@(xi)   ones(size(xi)) , valid_boundary_input);
            addOptional(p,'Gx4' ,@(eta) -ones(size(eta)), valid_boundary_input);
            addOptional(p,'Gy4' ,@(eta)              eta, valid_boundary_input);
            addOptional(p,'Gx1_domain',[-1,1],validVec);
            addOptional(p,'Gy1_domain',[-1,1],validVec);
            addOptional(p,'Gx2_domain',[-1,1],validVec);
            addOptional(p,'Gy2_domain',[-1,1],validVec);
            addOptional(p,'Gx3_domain',[-1,1],validVec);
            addOptional(p,'Gy3_domain',[-1,1],validVec);
            addOptional(p,'Gx4_domain',[-1,1],validVec);
            addOptional(p,'Gy4_domain',[-1,1],validVec);
            parse(p,varargin{:});
            this.G{1,1} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gx1, ...
                                                       p.Results.Gx1_domain(1),...
                                                       p.Results.Gx1_domain(2));
            this.G{1,2} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gy1, ...
                                                       p.Results.Gy1_domain(1),...
                                                       p.Results.Gy1_domain(2));
            this.G{2,1} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gx2, ...
                                                       p.Results.Gx2_domain(1),...
                                                       p.Results.Gx2_domain(2));
            this.G{2,2} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gy2, ...
                                                       p.Results.Gy2_domain(1),...
                                                       p.Results.Gy2_domain(2));
            this.G{3,1} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gx3, ...
                                                       p.Results.Gx3_domain(1),...
                                                       p.Results.Gx3_domain(2));
            this.G{3,2} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gy3, ...
                                                       p.Results.Gy3_domain(1),...
                                                       p.Results.Gy3_domain(2));
            this.G{4,1} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gx4, ...
                                                       p.Results.Gx4_domain(1),...
                                                       p.Results.Gx4_domain(2));
            this.G{4,2} = transfinite_interpolant_2D.setup_boundaries( ...
                                                       p.Results.Gy4, ...
                                                       p.Results.Gy4_domain(1),...
                                                       p.Results.Gy4_domain(2));
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
        function plot_grid(this,n,varargin)
            t = linspace(-1,1,n);
            [x,y] = this.eval_grid(t,t);
            hold on
            plot(x,y,'k'); plot(x.',y.','k');
            axis equal
        end
        function [x,y] = eval_grid(this,xi_vec,eta_vec)
            [xi,eta] = ndgrid(xi_vec,eta_vec);
            x = this.eval_x(xi,eta);
            y = this.eval_y(xi,eta);
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
        function F = setup_boundaries(G,t_min,t_max)
            if isa(G,'function_handle')
                F = @(xi) G( transfinite_interpolant_2D.inv_linear_map(xi,t_min,t_max) );
                return
            elseif isscalar(G)
                F = @(x) G*ones(size(x));
            else
                t = linspace(-1,1,numel(G));
                pp = csape(t,G);
                F = @(x) fnval(pp,x);
            end
        end
        function xi = linear_map(t,t_min,t_max)
            xi = 2*(t-t_min)/(t_max-t_min) - 1;
        end
        function t = inv_linear_map(xi,t_min,t_max)
            t = t_min + (t_max-t_min)*(xi+1)/2;
        end
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