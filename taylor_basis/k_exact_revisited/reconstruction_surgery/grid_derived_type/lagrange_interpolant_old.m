classdef lagrange_interpolant_old

    properties
        q (:,:)   {mustBeNumeric} = 0
        w (:,:)   {mustBeNumeric} = 0
        D (:,1)
    end

    methods
        function this = lagrange_interpolant_old(Nmax)
            if nargin < 1
                [this.q,this.w,this.D] = this.generate_1D_barycentric_info(10);
            else
                [this.q,this.w,this.D] = this.generate_1D_barycentric_info(Nmax);
            end
        end
        function J  = jacobian_determinant(this,X1,X2,X3,xq,dim)
            if dim == 1
                J = this.mapping_jacobian_determinant_1D(xq,X1,X2,X3,this.q,this.w,this.D);
            elseif dim == 2
                J = this.mapping_jacobian_determinant_2D(xq,X1,X2,X3,this.q,this.w,this.D);
            elseif dim == 3
                J = this.mapping_jacobian_determinant_3D(xq,X1,X2,X3,this.q,this.w,this.D);
            end
        end
        function Ja = metrics(this,X1,X2,X3,xq)
            Ja = this.mapping_curl_invariant_metrics(xq,X1,X2,X3,this.q,this.w,this.D);
        end
        function [a1,a2] = base_vectors(this,X1,X2,X3,xq)
            a1 = this.covariant_base_vectors_3D(xq,X1,X2,X3,this.q,this.w,this.D);
            a2 = this.contravariant_base_vectors_3D(xq,X1,X2,X3,this.q,this.w,this.D);
        end
        function kappa = normal_curvature(this,X1,X2,X3,xq,dir,vec)
            kappa = this.get_curvature(xq,dir,vec,X1,X2,X3,this.q,this.w,this.D);
        end
        function yq = map_point(this,X1,X2,X3,xq)
            N = size(X1,[1,2,3]);
            yq(1,1) = this.lagbary_3D(xq,this.q,X1,this.w,N);
            yq(2,1) = this.lagbary_3D(xq,this.q,X2,this.w,N);
            yq(3,1) = this.lagbary_3D(xq,this.q,X3,this.w,N);
        end
        function [Yq1,Yq2,Yq3] = map_grid(this,X1,X2,X3,N1,N2,N3)
            Yq1 = zeros(N1,N2,N3);
            Yq2 = zeros(N1,N2,N3);
            Yq3 = zeros(N1,N2,N3);
            if N3 > 1
                [Xq1,Xq2,Xq3] = ndgrid(linspace(-1,1,N1),linspace(-1,1,N2),linspace(-1,1,N3));
            else
                [Xq1,Xq2] = ndgrid(linspace(-1,1,N1),linspace(-1,1,N2));
                Xq3 = zeros(N1,N2);
            end
            
            for k = 1:N3
                for j = 1:N2
                    for i = 1:N1
                        xq = [Xq1(i,j,k),Xq2(i,j,k),Xq3(i,j,k)];
                        yq = this.map_point(X1,X2,X3,xq);
                        Yq1(i,j,k) = yq(1);
                        Yq2(i,j,k) = yq(2);
                        Yq3(i,j,k) = yq(3);
                    end
                end
            end
        end
    end

    methods (Static)
        function result = almost_equal(a,b)
            test1 = ( (a==0) | (b==0) );
            test2 = ( abs(a-b) <= 2*eps(1) );
            test3 = ( abs(a-b) <= eps(abs(a)) ) & ( abs(a-b) <= eps(abs(b)) );
            result = logical( test1.*test2 + (~test1).*test3 );
        end
        function w = barycentric_weights(x,N)
            w = ones(N,1);
            for j = 2:N       % for j = 1:N-1
                for k = 1:j-1 % for k = 0:j-1
                    w(k) = w(k) * ( x(k) - x(j) );
                    w(j) = w(j) * ( x(j) - x(k) );
                end
            end
            w = 1./w;
        end
        function D = polynomial_derivative_matrix(x,w,N)
            D = zeros(N,N);
            for i = 1:N
                for j = 1:N
                    if j ~= i
                        D(i,j) = w(j)/w(i) * 1 / (x(i) - x(j));
                        D(i,i) = D(i,i) - D(i,j);
                    end
                end
            end
        end
        function D = mth_order_polynomial_derivative_matrix(x,w,N,m)
            D = zeros(N,N,m);
            D(:,:,1) = lagrange_interpolant_old.polynomial_derivative_matrix(x,w,N);
            for k = 2:m
                for i = 1:N
                    D(i,i,k) = 0;
                    for j = 1:N
                        if j ~= i
                            D(i,j,k) = ( k / (x(i) - x(j)) )*(w(j)/w(i)*D(i,i,k-1) - D(i,j,k-1) );
                            D(i,i,k) = D(i,i,k) - D(i,j,k);
                        end
                    end
                end
            end
        end
        % function [q,w,D] = generate_1D_barycentric_info(Nmax)
        %     q = zeros(Nmax,Nmax);
        %     w = zeros(Nmax,Nmax);
        %     D = zeros(Nmax,Nmax,Nmax);
        %     q(1,1) = 0;
        %     w(1,1) = 1;
        %     D(1,1) = 0;
        %     for i = 2:Nmax
        %         q(1:i,i) = linspace(-1,1,i);
        %         w(1:i,i) = lagrange_interpolant_old.barycentric_weights(q(1:i,i),i);
        %         D(1:i,1:i,i) = lagrange_interpolant_old.polynomial_derivative_matrix( q(1:i,i), w(1:i,i), i );
        %     end
        % end
        function [q,w,D] = generate_1D_barycentric_info(Nmax)
            q = zeros(Nmax,Nmax);
            w = zeros(Nmax,Nmax);
            D = cell(Nmax,1);
            q(1,1) = 0;
            w(1,1) = 1;
            D{1}   = 0;
            for i = 2:Nmax
                q(1:i,i) = linspace(-1,1,i);
                w(1:i,i) = lagrange_interpolant_old.barycentric_weights(q(1:i,i),i);
                D{i} = lagrange_interpolant_old.mth_order_polynomial_derivative_matrix( q(1:i,i), w(1:i,i), i, 2 );
            end
        end
        function pval = lagbary(p,dir,f,q,w,Npts)
            % p - point in parametric space    1 x Ndim
            % dir - coordinate direction for interpolation
            % f - function values at quadrature points Npts
            % q - quadrature points                    Npts x Ndim
            % w - quadrature weights                   Npts x Ndim
            % Npts - number of quad pts in each direction
            A = 0; F = 0;
            % for each quadrature point in direction dir
            for i = 1:Npts(dir)
                x = q(i,Npts(dir)) - p;
                if lagrange_interpolant_old.almost_equal(x,0)
                    pval = f(i);
                    return;
                end
                t1 = w(i,Npts(dir))/x;
                A = A + t1 * f(i);
                F = F + t1;
            end
            pval = A/F;
        end
        function [pval,dpval] = lagbary_wderiv(p,dir,f,q,w,Dmat,Npts)
            % p - point in parametric space    1 x Ndim
            % dir - coordinate direction for derivative
            % f - function values at quadrature points Npts
            % q - quadrature points                    Npts x Ndim
            % w - quadrature weights                   Npts x Ndim
            % D - derivative matrices for quadrature pts
            % Npts - number of quad pts in each direction
            A = 0; B = 0; C = 0; F = 0;
            % for each quadrature point in direction dir
            for i = 1:Npts(dir)
                x = q(i,Npts(dir)) - p;
                if lagrange_interpolant_old.almost_equal(x,0)
                    pval = f(i);
                    % dpval = D(i,1:Npts(dir),Npts(dir)) * f(:);
                    dpval = Dmat{Npts(dir)}(i,1:Npts(dir),1) * f(:);
                    return;
                end
                t1 = w(i,Npts(dir))/x;
                A = A + t1 * f(i);
                F = F + t1;
                t2 = t1/x;
                B = B + t2 * f(i);
                C = C + t2;
            end
            pval = A/F;
            FF = F * F;
            AC = A * C;
            dpval = (B * F - AC) / FF;
        end
        function [val,dval,d2val] = lagbary_wderiv2(p,dir,f,q,w,Dmat,Npts)
            % p - point in parametric space    1 x Ndim
            % dir - coordinate direction for derivative
            % f - function values at quadrature points Npts
            % q - quadrature points                    Npts x Ndim
            % w - quadrature weights                   Npts x Ndim
            % D - derivative matrices for quadrature pts
            % Npts - number of quad pts in each direction
            A = 0; B = 0; C = 0; D = 0; E = 0; F = 0;
            % for each quadrature point in direction dir
            for i = 1:Npts(dir)
                x = q(i,Npts(dir)) - p;
                if lagrange_interpolant_old.almost_equal(x,0)
                    val   = f(i);
                    dval  = Dmat{Npts(dir)}(i,1:Npts(dir),1) * f(:);
                    d2val = Dmat{Npts(dir)}(i,1:Npts(dir),2) * f(:);
                    return;
                end
                t1 = w(i,Npts(dir))/x;
                A = A + t1 * f(i);
                F = F + t1;

                t2 = t1/x;
                B = B + t2 * f(i);
                C = C + t2;

                t3 = t2/x;
                D = D + t3 * f(i);
                E = E + t3;
            end
            val = A/F;

            FF = F * F;
            AC = A * C;
            dval  =  (B * F - AC) / FF;
            d2val =     ( 2 * D ) / F ...
                  - ( 2 * E * A ) / FF ...
                   - ( 2* B * C ) / FF ...
                 + ( 2 * C * AC ) / ( FF * F );
        end
        function val = lagbary_2D(p,q,f,w,N)
            tmp   = zeros(N(2),1);
            for j = 1:N(2)
                tmp(j) = lagrange_interpolant_old.lagbary(p(1),1,f(:,j),q,w,N);
            end
            val = lagrange_interpolant_old.lagbary(p(2),2,tmp,q,w,N);
        end
        function [val,grad] = lagbary_2D_wgrad(p,q,f,w,D,N)
            grad = zeros(2,1);
            tmp  = zeros(N(2),1);
            gtmp = zeros(N(2),1);
            for j = 1:N(2)
                [tmp(j),gtmp(j)] = lagrange_interpolant_old.lagbary_wderiv(p(1),1,f(:,j),q,w,D,N);
            end
            [val,grad(2)] = lagrange_interpolant_old.lagbary_wderiv(p(2),2,tmp,q,w,D,N);
            grad(1)   = lagrange_interpolant_old.lagbary(p(2),2,gtmp,q,w,N);
        end
        function [val,grad,hess] = lagbary_2D_whess(p,q,f,w,D,N)
            grad = zeros(2,1);
            hess = zeros(2,2);
            tmp  = zeros(N(2),1);
            gtmp = zeros(N(2),1);
            htmp = zeros(N(2),1);
            for j = 1:N(2)
                [ tmp(j), gtmp(j), htmp(j) ] = lagrange_interpolant_old.lagbary_wderiv2( p(1), 1, f(:,j),q,w,D,N);
            end
            [ val, grad(2), hess(2,2) ] = lagrange_interpolant_old.lagbary_wderiv2( p(2),2, tmp,q,w,D,N);
            [      grad(1), hess(1,2) ] = lagrange_interpolant_old.lagbary_wderiv(  p(2),2,gtmp,q,w,D,N);
            hess(1,1)                   = lagrange_interpolant_old.lagbary(         p(2),2,htmp,q,w,N);
            hess(2,1) = hess(1,2);
        end
        function val = lagbary_3D(p,q,f,w,N)
            tmp   = zeros(N(2),N(3));
            for k = 1:N(3)
                for j = 1:N(2)
                    tmp(j,k) = lagrange_interpolant_old.lagbary(p(1),1,f(:,j,k),q,w,N);
                end
            end
            for k = 1:N(3)
                tmp(1,k) = lagrange_interpolant_old.lagbary(p(2),2,tmp(:,k),q,w,N);
            end
            val = lagrange_interpolant_old.lagbary(p(3),3,tmp(1,:),q,w,N);
        end
        function [val,grad] = lagbary_3D_wgrad(p,q,f,w,D,N)
            grad = zeros(3,1);
            tmp   = zeros(N(2),N(3));
            gtmp1 = zeros(N(2),N(3));
            for k = 1:N(3)
                for j = 1:N(2)
                    [tmp(j,k),gtmp1(j,k)] = lagrange_interpolant_old.lagbary_wderiv(p(1),1,f(:,j,k),q,w,D,N);
                end
            end
            gtmp2 = zeros(1,N(3));
            for k = 1:N(3)
                [tmp(1,k),gtmp2(1,k)] = lagrange_interpolant_old.lagbary_wderiv(p(2),2,tmp(:,k),q,w,D,N);
                gtmp1(1,k)            = lagrange_interpolant_old.lagbary(p(2),2,gtmp1(:,k),q,w,N);
            end
            [val,grad(3)] = lagrange_interpolant_old.lagbary_wderiv(p(3),3,tmp(1,:),q,w,D,N);
            grad(2)       = lagrange_interpolant_old.lagbary(p(3),3,gtmp2(1,:),q,w,N);
            grad(1)       = lagrange_interpolant_old.lagbary(p(3),3,gtmp1(1,:),q,w,N);
        end
        function [val,grad,hess] = lagbary_3D_whess(p,q,f,w,D,N)
            grad = zeros(3,1);
            hess = zeros(3,3);
            tmp   = zeros(N(2),N(3));
            gtmp1 = zeros(N(2),N(3));
            htmp1 = zeros(N(2),N(3));
            for k = 1:N(3)
                for j = 1:N(2)
                    [ tmp(j,k), gtmp1(j,k), htmp1(j,k) ] = lagrange_interpolant_old.lagbary_wderiv2( p(1), 1, f(:,j,k),q,w,D,N);
                end
            end
            gtmp2 = zeros(1,N(3));
            htmp2 = zeros(1,N(3));
            htmp3 = zeros(1,N(3));
            for k = 1:N(3)
                [ tmp(1,k), gtmp2(1,k), htmp3(1,k)] = lagrange_interpolant_old.lagbary_wderiv2( p(2), 2,   tmp(:,k),q,w,D,N);
                [           gtmp1(1,k), htmp2(1,k)] = lagrange_interpolant_old.lagbary_wderiv(  p(2), 2, gtmp1(:,k),q,w,D,N);
                htmp1(1,k)                          = lagrange_interpolant_old.lagbary(         p(2), 2, htmp1(:,k),q,w,N);
            end
            [ val, grad(3), hess(3,3) ] = lagrange_interpolant_old.lagbary_wderiv2( p(3), 3,   tmp(1,:),q,w,D,N);
            [      grad(1), hess(1,3) ] = lagrange_interpolant_old.lagbary_wderiv(  p(3), 3, gtmp1(1,:),q,w,D,N);
            [      grad(2), hess(2,3) ] = lagrange_interpolant_old.lagbary_wderiv(  p(3), 3, gtmp2(1,:),q,w,D,N);
            hess(2,2)                   = lagrange_interpolant_old.lagbary(         p(3), 3, htmp3(1,:),q,w,N);
            hess(1,1)                   = lagrange_interpolant_old.lagbary(         p(3), 3, htmp1(1,:),q,w,N);
            hess(1,2)                   = lagrange_interpolant_old.lagbary(         p(3), 3, htmp2(1,:),q,w,N);

            hess(2,1) = hess(1,2);
            hess(3,1) = hess(1,3);
            hess(3,2) = hess(2,3);
        end

        function a = covariant_base_vectors_1D(point,X1,X2,X3,qn,wn,Dn)
            % tangent vector
            % X1, X2, X3 are 1D arrays
            Nn = size(X1,1);
            [~,a(1)] = lagrange_interpolant_old.lagbary_wderiv(point(1),1,X1,qn,wn,Dn,Nn);
            [~,a(2)] = lagrange_interpolant_old.lagbary_wderiv(point(2),1,X2,qn,wn,Dn,Nn);
            [~,a(3)] = lagrange_interpolant_old.lagbary_wderiv(point(3),1,X3,qn,wn,Dn,Nn);
        end
        function [t,n,b,kappa] = frenet_frame(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 1D arrays
            Nn = size(X1,1);

            % get tangent vector and derivative at point
            [~,t(1),tdot(1)] = lagrange_interpolant_old.lagbary_wderiv2(point(1),1,X1,qn,wn,Dn,Nn);
            [~,t(2),tdot(2)] = lagrange_interpolant_old.lagbary_wderiv2(point(2),1,X2,qn,wn,Dn,Nn);
            [~,t(3),tdot(3)] = lagrange_interpolant_old.lagbary_wderiv2(point(3),1,X3,qn,wn,Dn,Nn);

            % calculate kappa
            tmp = cross(t,tdot);
            len = sqrt(sum(t.^2));
            kappa = sqrt(sum(tmp.^2))/len^3;

            % get unit tangent
            t = t/sqrt(sum(t.^2));
            
            tol = 1e-12;
            if kappa>tol
                % get unit tangent along edge
                t1 = zeros(Nn,3);
                for i = 1:Nn
                    [~,t1(i,1)] = lagrange_interpolant_old.lagbary_wderiv(X1(i),1,X1,qn,wn,Dn,Nn);
                    [~,t1(i,2)] = lagrange_interpolant_old.lagbary_wderiv(X2(i),1,X2,qn,wn,Dn,Nn);
                    [~,t1(i,3)] = lagrange_interpolant_old.lagbary_wderiv(X3(i),1,X3,qn,wn,Dn,Nn);
                    t1(i,:) = t1(i,:)/sqrt(sum(t1(i,:).^2));
                end
                % get derivative of unit tangent at point
                [~,t1dot(1)] = lagrange_interpolant_old.lagbary_wderiv(point(1),1,t1(:,1),qn,wn,Dn,Nn);
                [~,t1dot(2)] = lagrange_interpolant_old.lagbary_wderiv(point(2),1,t1(:,2),qn,wn,Dn,Nn);
                [~,t1dot(3)] = lagrange_interpolant_old.lagbary_wderiv(point(3),1,t1(:,3),qn,wn,Dn,Nn);
                n = t1dot/sqrt(sum(t1dot.^2));
                b = cross(t,n);
            else
                % find direction with minimum variance
                [~,i] = min([sum(X1.^2),sum(X2.^2),sum(X3.^2)]);
                b = zeros(3,1);
                b(i) = 1;
                n = cross(b,t);
            end
        end
        function a = covariant_metric_tensor_1D(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 1D arrays
            a1 = lagrange_interpolant_old.covariant_base_vectors_1D(point,X1,X2,X3,qn,wn,Dn);
            a = dot(a1,a1);
        end
        function a = contravariant_metric_tensor_1D(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 1D arrays
            a1 = lagrange_interpolant_old.covariant_metric_tensor_1D(point,X1,X2,X3,qn,wn,Dn);
            a = 1/a1;
        end
        function J = mapping_jacobian_determinant_1D(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 1D arrays
            a = lagrange_interpolant_old.covariant_metric_tensor_1D(point,X1,X2,X3,qn,wn,Dn);
            J = sqrt(a);
        end
        function a = covariant_base_vectors_2D(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 2D arrays
            Nn = size(X1,[1,2]);
            a = zeros(3,2);
            [~,a(1,:)] = lagrange_interpolant_old.lagbary_2D_wgrad(point,qn,X1,wn,Dn,Nn);
            [~,a(2,:)] = lagrange_interpolant_old.lagbary_2D_wgrad(point,qn,X2,wn,Dn,Nn);
            [~,a(3,:)] = lagrange_interpolant_old.lagbary_2D_wgrad(point,qn,X3,wn,Dn,Nn);
        end
        function a = covariant_metric_tensor_2D(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 2D arrays
            a1 = lagrange_interpolant_old.covariant_base_vectors_2D(point,X1,X2,X3,qn,wn,Dn);
            a = zeros(2,2);
            a(1,1) = dot(a1(:,1),a1(:,1));
            a(1,2) = dot(a1(:,1),a1(:,2));
            a(2,1) = a(1,2);
            a(2,2) = dot(a1(:,2),a1(:,2));
        end
        function a = contravariant_metric_tensor_2D(point,X1,X2,X3,qn,wn,Dn)
            % X1, X2, X3 are 2D arrays
            a1 = covariant_metric_tensor_2D(point,X1,X2,X3,qn,wn,Dn);
            sa = sum(cross(a1(:,1),a1(:,2)).^2);
            a = zeros(2,2);
            a(1,1) =  a1(2,2)/sa;
            a(1,2) = -a1(1,2)/sa;
            a(2,1) =  a(2,1);
            a(2,2) =  a1(1,1)/sa;
        end

        function a = covariant_base_vectors_3D(point,X1,X2,X3,qn,wn,Dn)
            Nn = size(X1,[1,2,3]);
            a = zeros(3,3);
            [~,a(1,:)] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X1,wn,Dn,Nn);
            [~,a(2,:)] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X2,wn,Dn,Nn);
            [~,a(3,:)] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X3,wn,Dn,Nn);
        end

        function a = contravariant_base_vectors_2D(point,X1,X2,X3,qn,wn,Dn)
            Ja = lagrange_interpolant_old.mapping_curl_invariant_metrics(point,X1,X2,X3,qn,wn,Dn);
            J  = lagrange_interpolant_old.mapping_jacobian_determinant_3D(point,X1,X2,X3,qn,wn,Dn);
            a = Ja/J;
        end
        function a = contravariant_base_vectors_3D(point,X1,X2,X3,qn,wn,Dn)
            Ja = lagrange_interpolant_old.mapping_curl_invariant_metrics(point,X1,X2,X3,qn,wn,Dn);
            J  = lagrange_interpolant_old.mapping_jacobian_determinant_3D(point,X1,X2,X3,qn,wn,Dn);
            a = Ja/J;
        end
        function J = mapping_jacobian_determinant_2D(point,X1,X2,~,qn,wn,Dn)
            Nn = size(X1,[1,2,3]);
            Nn(3) = 1;
            X1 = X1(:,:,1);
            X2 = X2(:,:,1);
            a = zeros(2,2);
            [~,tmp] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X1,wn,Dn,Nn);
            a(1,:) = tmp(1:2);
            [~,tmp] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X2,wn,Dn,Nn);
            a(2,:) = tmp(1:2);
            J = det(a);
        end
        function J = mapping_jacobian_determinant_3D(point,X1,X2,X3,qn,wn,Dn)
            Nn = size(X1,[1,2,3]);
            a = zeros(3,3);
            [~,a(1,:)] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X1,wn,Dn,Nn);
            [~,a(2,:)] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X2,wn,Dn,Nn);
            [~,a(3,:)] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X3,wn,Dn,Nn);
            J = det(a);
        end
        function Ja = mapping_curl_invariant_metrics(point,X1,X2,X3,qn,wn,Dn)
            Nn = size(X1,[1,2,3]);
            X = cat(4,X1,X2,X3);
            Ja = zeros(3,3); % basis (i) (1-3), basis component (n) (1-3)
            nml = 1:3;
            for n = 1:3
                nml2 = mod(nml+n+1,3)+1;
                X_l = X(:,:,:,nml2(3)); % Nn(1) x Nn(2) x Nn(3)
                X_m = X(:,:,:,nml2(2)); % Nn(1) x Nn(2) x Nn(3)
                [~,dX_l] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X_l,wn,Dn,Nn); % 3 x 1
                [~,dX_m] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,X_m,wn,Dn,Nn); % 3 x 1
                d1 = X_l*dX_m(1) - X_m*dX_l(1); % Nn(1) x Nn(2) x Nn(3)
                d2 = X_l*dX_m(2) - X_m*dX_l(2); % Nn(1) x Nn(2) x Nn(3)
                d3 = X_l*dX_m(3) - X_m*dX_l(3); % Nn(1) x Nn(2) x Nn(3)
                [~,dd1] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,d1,wn,Dn,Nn); % 3 x 1
                [~,dd2] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,d2,wn,Dn,Nn); % 3 x 1
                [~,dd3] = lagrange_interpolant_old.lagbary_3D_wgrad(point,qn,d3,wn,Dn,Nn); % 3 x 1
                Ja(n,1) = -0.5*( dd3(2) - dd2(3) );
                Ja(n,2) = -0.5*( dd1(3) - dd3(1) );
                Ja(n,3) = -0.5*( dd2(1) - dd1(2) );
            end
        end

        function kappa = get_curvature(point,dir,vec,X1,X2,X3,qn,wn,Dn)
            vec = vec/sqrt(dot(vec,vec));
            Nn = size(X1,[1,2,3]);
            a1 = lagrange_interpolant_old.contravariant_base_vectors_3D(point,X1,X2,X3,qn,wn,Dn);
            hess = zeros(3,3,3);
            [~,~,hess(:,:,1)] = lagrange_interpolant_old.lagbary_3D_whess(point,qn,X1,wn,Dn,Nn);
            [~,~,hess(:,:,2)] = lagrange_interpolant_old.lagbary_3D_whess(point,qn,X2,wn,Dn,Nn);
            [~,~,hess(:,:,3)] = lagrange_interpolant_old.lagbary_3D_whess(point,qn,X3,wn,Dn,Nn);
            if dir == 1
                n = a1(:,1)/sqrt(sum(a1(:,1).^2));
                vec = vec*( 1 - dot(vec,n)/sqrt(dot(vec,vec)) );
                kappa = dot( n, squeeze( hess(2,2,:) ) )*dot(vec,a1(:,2))^2 ...
                    + 2*dot( n, squeeze( hess(2,3,:) ) )*dot(vec,a1(:,2))*dot(vec,a1(:,3)) ...
                      + dot( n, squeeze( hess(3,3,:) ) )*dot(vec,a1(:,3))^2;
            elseif dir == 2
                n = a1(:,2)/sqrt(sum(a1(:,2).^2));
                vec = vec*( 1 - dot(vec,n)/sqrt(dot(vec,vec)) );
                kappa = dot( n, squeeze( hess(1,1,:) ) )*dot(vec,a1(:,1))^2 ...
                    + 2*dot( n, squeeze( hess(1,3,:) ) )*dot(vec,a1(:,1))*dot(vec,a1(:,3)) ...
                      + dot( n, squeeze( hess(3,3,:) ) )*dot(vec,a1(:,3))^2;
            else
                n = a1(:,3)/sqrt(sum(a1(:,3).^2));
                vec = vec*( 1 - dot(vec,n)/sqrt(dot(vec,vec)) );
                kappa = dot( n, squeeze( hess(1,1,:) ) )*dot(vec,a1(:,1))^2 ...
                    + 2*dot( n, squeeze( hess(1,2,:) ) )*dot(vec,a1(:,1))*dot(vec,a1(:,2)) ...
                      + dot( n, squeeze( hess(2,2,:) ) )*dot(vec,a1(:,2))^2;
            end
        end
    end

end