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
            validRealPoints  = @(x) all(arrayfun(@(x)validReal(x),x),"all") && isvector(x);
            validVec         = @(x) validRealPoints(x) && numel(x)==2;
            valid_fcn_handle = @(x) isa(x,'function_handle');
            valid_boundary_input = @(x) valid_fcn_handle(x) || validRealPoints(x);
            addOptional(p,'G1x' ,@(xi)                xi, valid_boundary_input);
            addOptional(p,'G1y' ,@(xi)  -ones(size(xi)) , valid_boundary_input);
            addOptional(p,'G2x' ,@(eta)  ones(size(eta)), valid_boundary_input);
            addOptional(p,'G2y' ,@(eta)              eta, valid_boundary_input);
            addOptional(p,'G3x' ,@(xi)                xi, valid_boundary_input);
            addOptional(p,'G3y' ,@(xi)   ones(size(xi)) , valid_boundary_input);
            addOptional(p,'G4x' ,@(eta) -ones(size(eta)), valid_boundary_input);
            addOptional(p,'G4y' ,@(eta)              eta, valid_boundary_input);
            addOptional(p,'G1x_dom',[-1,1],validVec);
            addOptional(p,'G1y_dom',[-1,1],validVec);
            addOptional(p,'G2x_dom',[-1,1],validVec);
            addOptional(p,'G2y_dom',[-1,1],validVec);
            addOptional(p,'G3x_dom',[-1,1],validVec);
            addOptional(p,'G3y_dom',[-1,1],validVec);
            addOptional(p,'G4x_dom',[-1,1],validVec);
            addOptional(p,'G4y_dom',[-1,1],validVec);
            parse(p,varargin{:});
            results = transfinite_interpolant_2D.check_boundaries(p);
            this.G{1,1} = results.G1x;
            this.G{1,2} = results.G1y;
            this.G{2,1} = results.G2x;
            this.G{2,2} = results.G2y;
            this.G{3,1} = results.G3x;
            this.G{3,2} = results.G3y;
            this.G{4,1} = results.G4x;
            this.G{4,2} = results.G4y;
            % this.G = cellfun(@(x)transfinite_interpolant_2D.setup_boundaries(x),this.G,'UniformOutput',false);
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
        function results = check_boundaries(p)
            results = structfun(@(x)transfinite_interpolant_2D.setup_boundaries(x),p.Results,'UniformOutput',false);
            % continuity
            % G4( 1) == G1(-1)
            % G1( 1) == G2(-1)
            % G2( 1) == G3(-1)
            % G3( 1) == G4(-1)
            % xy_eval = zeros(8,2);
            % tol = 1.0e-12;
            % 
            % xy_eval(1,1) = results.G1x(-1); xy_eval(1,2) = results.G4x( 1);
            % xy_eval(2,1) = results.G1y(-1); xy_eval(2,2) = results.G4y( 1);
            % xy_eval(3,1) = results.G2x(-1); xy_eval(3,2) = results.G1x( 1);
            % xy_eval(4,1) = results.G2y(-1); xy_eval(4,2) = results.G1y( 1);
            % xy_eval(5,1) = results.G3x(-1); xy_eval(5,2) = results.G2x( 1);
            % xy_eval(6,1) = results.G3y(-1); xy_eval(6,2) = results.G2y( 1);
            % xy_eval(7,1) = results.G4x(-1); xy_eval(7,2) = results.G3x( 1);
            % xy_eval(8,1) = results.G4y(-1); xy_eval(8,2) = results.G3y( 1);
            % 
            % % if it all matches up, no need to change anything
            % matches = abs(xy_eval(:,1)-xy_eval(:,2))<tol;
            % if all(matches)
            %     return
            % end
            % 
            % % otherwise
            % 
            % % find which ones are user input
            % n = numel(p.Parameters);
            % user_input = false(n,1);
            % for i = 1:n
            %     if ( ~any(cellfun(@(s2)strcmp(p.Parameters{i},s2),p.UsingDefaults)) )
            %         user_input(i) = true;
            %     end
            % end
            % m = sum(user_input);
            % 
            % % this.G{1,1} = results.G1x;
            % % this.G{1,2} = results.G1y;
            % % this.G{2,1} = results.G2x;
            % % this.G{2,2} = results.G2y;
            % % this.G{3,1} = results.G3x;
            % % this.G{3,2} = results.G3y;
            % % this.G{4,1} = results.G4x;
            % % this.G{4,2} = results.G4y;
            % 
            % flags = false(n,2); % lo or hi mismatch
            % vals  = [-ones(n,1);ones(n,1)]; % lo or hi value
            % % continuity
            % % G4( 1) == G1(-1)
            % % G1( 1) == G2(-1)
            % % G2( 1) == G3(-1)
            % % G3( 1) == G4(-1)
            % % G1x(-1)
            % if ~matches(1)
            %     if user_input(1) && user_input(7)
            %         error('G1x(-1) must equal G4x(1)')
            %     elseif ~user_input(1)
            %         flags(1,1) = true;
            %     else
            %         flags(7,2) = true;
            %     end
            % end
            % 
            % % G1y(-1)
            % if ~matches(2)
            %     if user_input(2) && user_input(8)
            %         error('G1y(-1) must equal G4y(1)')
            %     elseif ~user_input(2)
            %         flags(2,1) = true;
            %     else
            %         flags(8,2) = true;
            %     end
            % end
            % 
            % % G2x(-1)
            % if ~matches(3)
            %     if user_input(3) && user_input(1)
            %         error('G2x(-1) must equal G1x(1)')
            %     elseif ~user_input(3)
            %         flags(3,1) = true;
            %     else
            %         flags(1,2) = true;
            %     end
            % end
            % 
            % % G2y(-1)
            % if ~matches(4)
            %     if user_input(4) && user_input(2)
            %         error('G2y(-1) must equal G1y(1)')
            %     elseif ~user_input(4)
            %         flags(4,1) = true;
            %     else
            %         flags(2,2) = true;
            %     end
            % end
            % 
            % % G3x(-1)
            % if ~matches(5)
            %     if user_input(5) && user_input(3)
            %         error('G3x(-1) must equal G2x(1)')
            %     elseif ~user_input(5)
            %         flags(5,1) = true;
            %     else
            %         flags(3,2) = true;
            %     end
            % end
            % 
            % % G3y(-1)
            % if ~matches(6)
            %     if user_input(6) && user_input(4)
            %         error('G3y(-1) must equal G2y(1)')
            %     elseif ~user_input(6)
            %         flags(6,1) = true;
            %     else
            %         flags(4,2) = true;
            %     end
            % end
            % 
            % % G4x(-1)
            % if ~matches(7)
            %     if user_input(7) && user_input(5)
            %         error('G4x(-1) must equal G3x(1)')
            %     elseif ~user_input(7)
            %         flags(7,1) = true;
            %     else
            %         flags(5,2) = true;
            %     end
            % end
            % 
            % % G4y(-1)
            % if ~matches(8)
            %     if user_input(8) && user_input(6)
            %         error('G4y(-1) must equal G3y(1)')
            %     elseif ~user_input(8)
            %         flags(8,1) = true;
            %     else
            %         flags(6,2) = true;
            %     end
            % end
            % 
            % for i = 1:n
            % 
            %     if any(flags(i,:))
            %         flaglo = flags(i,1);
            %         flaghi = flags(i,2);
            %         lo = -1*(~flaglo) + xy_eval(i,2);
            %             hi =  1;
            %             results.(p.Parameters{i}) = 

            % if ()
            % if abs(p.Results.G1x(-1)-p.Results.G1x(-1)
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