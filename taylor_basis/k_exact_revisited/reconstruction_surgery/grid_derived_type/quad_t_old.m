classdef quad_t_old
    properties
        n_quad    (1,1) {mustBeInteger} = 0  % # of quadrature points
        quad_pts  (:,:) {mustBeNumeric}      % Quadrature points in physical space
        quad_wts  (:,1) {mustBeNumeric}      % Quadrature weights
    end

    methods
        function this = allocate_quad( this, n_quad )
            this.n_quad = n_quad;
            this.quad_pts = zeros(3,n_quad);
            this.quad_wts = zeros(n_quad,1);
        end

        function integral = integrate(this,f)
            integral = 0;
            for n = 1:this.n_quad
                integral = integral + f(:,n)*this.quad_wts(n);
            end
        end

    end % methods

    methods (Static)
        function quad_ref = create_quad_ref_1D(n_quad)
            quad_ref = quad_t_old();
            quad_ref = quad_ref.allocate_quad( n_quad );
            [pts_1D,wts_1D] = quad_t_old.LegendreGaussNodesAndWeights(n_quad);
            quad_ref.quad_pts(1,:) = pts_1D;
            quad_ref.quad_wts      = wts_1D;
        end
        function quad_ref = create_quad_ref_2D(n_quad)
            quad_ref = quad_t_old();
            quad_ref = quad_ref.allocate_quad( n_quad^2 );
            [pts_1D,wts_1D] = quad_t_old.LegendreGaussNodesAndWeights(n_quad);
            cnt = 0;
            for j = 1:n_quad
                for i = 1:n_quad
                    cnt = cnt + 1;
                    quad_ref.quad_pts(:,cnt) = [pts_1D(i),pts_1D(j),0];
                    quad_ref.quad_wts(cnt)   = wts_1D(i)*wts_1D(j);
                end
            end
        end
        function quad_ref = create_quad_ref_3D(n_quad)
            quad_ref = quad_t_old();
            quad_ref = quad_ref.allocate_quad( n_quad^3 );
            [pts_1D,wts_1D] = quad_t_old.LegendreGaussNodesAndWeights(n_quad);
            cnt = 0;
            for k = 1:n_quad
                for j = 1:n_quad
                    for i = 1:n_quad
                        cnt = cnt + 1;
                        quad_ref.quad_pts(:,cnt) = [pts_1D(i),pts_1D(j),pts_1D(k)];
                        quad_ref.quad_wts(cnt)   = wts_1D(i)*wts_1D(j)*wts_1D(k);
                    end
                end
            end
        end
        function quad_phys = map_quad_ref_1D(X1,X2,X3,quad_ref,interpolant)
            quad_phys = quad_t_old();
            quad_phys = quad_phys.allocate_quad(quad_ref.n_quad);
            for n = 1:quad_ref.n_quad
                quad_phys.quad_pts(:,n) = interpolant.map_point(X1,X2,X3,quad_ref.quad_pts(:,n));
                det_jac = interpolant.jacobian_determinant(X1,X2,X3,quad_ref.quad_pts(:,n),1);
                quad_phys.quad_wts(n) = det_jac * quad_ref.quad_wts(n);
            end
        end
        function quad_phys = map_quad_ref_2D(X1,X2,X3,quad_ref,interpolant)
            quad_phys = quad_t_old();
            quad_phys = quad_phys.allocate_quad(quad_ref.n_quad);
            for n = 1:quad_ref.n_quad
                quad_phys.quad_pts(:,n) = interpolant.map_point(X1,X2,X3,quad_ref.quad_pts(:,n));
                det_jac = interpolant.jacobian_determinant(X1,X2,X3,quad_ref.quad_pts(:,n),2);
                quad_phys.quad_wts(n) = det_jac * quad_ref.quad_wts(n);
            end
        end
        function quad_phys = map_quad_ref_3D(X1,X2,X3,quad_ref,interpolant)
            quad_phys = quad_t_old();
            quad_phys = quad_phys.allocate_quad(quad_ref.n_quad);
            for n = 1:quad_ref.n_quad
                quad_phys.quad_pts(:,n) = interpolant.map_point(X1,X2,X3,quad_ref.quad_pts(:,n));
                det_jac = interpolant.jacobian_determinant(X1,X2,X3,quad_ref.quad_pts(:,n),3);
                quad_phys.quad_wts(n) = det_jac * quad_ref.quad_wts(n);
            end
        end
        function [x,w] = LegendreGaussNodesAndWeights(N)
            N = N-1;
            x = zeros(N+1,1);
            w = zeros(N+1,1);
            tol  = 4*eps(1);
            n_it = 10;
            if N == 0
                x(1) = 0;
                w(1) = 2;
            elseif N == 1
                x(1) = -sqrt(1/3);
                w(1) = 1;
                x(2) = -x(1);
                w(2) = w(1);
            else
                for j = 0:floor((N+1)/2)-1
                    x(j+1) = -cos( ( (2*j+1)/(2*N+2) )*pi );
                    for k = 1:n_it
                        [LNp1,dLNp1] = quad_t_old.LegendrePolynomialAndDerivative(N+1,x(j+1));
                        delta = -LNp1/dLNp1;
                        x(j+1) = x(j+1) + delta;
                        if abs(delta)<=tol*abs(x(j+1))
                            break;
                        end
                    end
                    [~,dLNp1] = quad_t_old.LegendrePolynomialAndDerivative(N+1,x(j+1));
                    x(N+1-j) = -x(j+1);
                    w(j+1) = 2/( (1-x(j+1)^2)*dLNp1^2);
                    w(N+1-j) = w(j+1);
                end
                if mod(N,2)==0
                    [~,dLNp1] = quad_t_old.LegendrePolynomialAndDerivative(N+1,0.0);
                    x(N/2+1) = 0;
                    w(N/2+1) = 2/dLNp1^2;
                end
            end
        end

        function [LN,dLN] = LegendrePolynomialAndDerivative(N,x)
            if N == 0
                LN  = 1;
                dLN = 0;
            elseif N == 1
                LN  = x;
                dLN = 1;
            else
                LNm2 = 1;
                LNm1 = x;
                dLNm2 = 0;
                dLNm1 = 1;
                for k = 2:N
                    LN = ( (2*k-1)/k )*x*LNm1 - ( (k-1)/k )*LNm2;
                    dLN   = dLNm2 + (2*k-1)*LNm1;
                    LNm2  = LNm1;
                    LNm1  = LN;
                    dLNm2 = dLNm1;
                    dLNm1 = dLN;
                end
            end
        end

    end

end