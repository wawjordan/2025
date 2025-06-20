classdef reference_quad_old
    properties
        Nq (3,1) {mustBeNumeric}
        e  (:,1) sub_quad_old
        f  (:,1) sub_quad_old
        v  (:,1) sub_quad_old
    end
    methods
        function this = reference_quad_old(Nq)
            if nargin > 0
                if numel(Nq)==2
                    this = this.create_reference_quad_2D(Nq);
                else
                    this = this.create_reference_quad_3D(Nq);
                end
            end
        end
        function this = create_reference_quad_2D(this,Nq)
            this.e = repmat(sub_quad_old(),0,1);
            this.f = repmat(sub_quad_old(),4,1);
            this.v = repmat(sub_quad_old(),1,1);
            this.Nq(1) = Nq(1);
            this.Nq(2) = Nq(2);
            this.Nq(3) = 1;
            [x1D,w1D] = this.generate_1D_quads(max(Nq));
            for i = 1:4
                [this.f(i).pts,this.f(i).wts] = this.edge_quad_2D( this.Nq, i, x1D, w1D );
            end
            [this.v(1).pts,this.v(1).wts] = this.face_quad_2D( this.Nq, 1, x1D, w1D );
        end
        function this = create_reference_quad_3D(this,Nq)
            this.e = repmat(sub_quad_old(),12,1);
            this.f = repmat(sub_quad_old(),6,1);
            this.v = repmat(sub_quad_old(),1,1);
            this.Nq(1) = Nq(1);
            this.Nq(2) = Nq(2);
            this.Nq(3) = Nq(3);
            [x1D,w1D] = this.generate_1D_quads(max(Nq));

            for i = 1:12
                [this.e(i).pts,this.e(i).wts] = this.edge_quad_3D( this.Nq, i, x1D, w1D );
            end
            for i = 1:6
                [this.f(i).pts,this.f(i).wts] = this.face_quad_3D( this.Nq, i, x1D, w1D );
            end
            [this.f(1).pts,this.f(1).wts] = this.volm_quad_3D( this.Nq, x1D, w1D );
        end
    end

    methods (Static)
        function [x,w] = generate_1D_quads(Nmax)
            % Nmax - maximum number of quadrature points in 1 dim
            % x    - quadrature abscissae
            % w    - quadrature weights
            x = zeros(Nmax,Nmax);
            w = zeros(Nmax,Nmax);
            for i = 1:Nmax
                [x(1:i,i),w(1:i,i)] = quad_t_old.LegendreGaussNodesAndWeights(i);
            end
        end
        function b = get_end_pts(i,dim)
            b = zeros(1,dim);
            tmp = i-1;
            for k = 1:dim
                b(k) = mod(tmp,2);
                tmp = fix(tmp/2);
            end
        end
        function [pts,wts] = edge_quad_2D(N,edge,x1D,w1D)
            pts = zeros(2,N(fix((edge-1)/2)+1));
            dir1 = fix( (edge-1)/2 ) + 1;
            dir2 = mod(dir1,2) + 1;
            Npts = N(dir1);
            end_pt = reference_quad_old.get_end_pts(mod(edge-1,2)+1,1);
            pts(dir1,:) = x1D(1:Npts,N(dir1));
            pts(dir2,:) = end_pt(1);
            wts         = w1D(1:Npts,N(dir1));
        end
        function [pts,wts] = edge_quad_3D(N,edge,x1D,w1D)
            pts = zeros( 3, N( fix( (edge-1)/4 ) + 1 ) );
            dir1 = fix( (edge-1)/4 ) + 1;
            dir2 = mod(dir1,3) + 1;
            dir3 = mod(dir1+1,3) + 1;
            Npts = N(dir1);
            if (dir2 > dir3)
                tmp = dir3;
                dir3 = dir2;
                dir2 = tmp;
            end
            end_pts = reference_quad_old.get_end_pts(mod(edge-1,4)+1,2);
            pts(dir1,:) = x1D(1:Npts,N(dir1));
            pts(dir2,:) = end_pts(1);
            pts(dir3,:) = end_pts(2);
            wts         = w1D(1:Npts,N(dir1));
        end
        function [pts,wts] = face_quad_2D(N,~,x1D,w1D)
            pts = zeros( 2, prod(N) );
            wts = zeros( prod(N), 1 );
            m = 0;
            for j = 1:N(2)
                for i = 1:N(1)
                    m = m + 1;
                    pts(1,m) = x1D(i,N(1));
                    pts(2,m) = x1D(j,N(2));
                    wts(m)   = w1D(i,N(1))*w1D(j,N(2));
                end
            end
        end
        function [pts,wts] = face_quad_3D(N,face,x1D,w1D)
            dir1 = fix( (face-1)/2 ) + 1;
            dir2 = mod(dir1  ,3) + 1;
            dir3 = mod(dir1+1,3) + 1;
            pts = zeros( 3, N(dir2)*N(dir3) );
            wts = zeros( N(dir2)*N(dir3), 1 );
            if (dir2 > dir3)
                tmp = dir3;
                dir3 = dir2;
                dir2 = tmp;
            end
            end_pt = reference_quad_old.get_end_pts(mod(face-1,2)+1,1);
            end_pt(1) = end_pt(1) * ( N(dir1) > 1 );
            m = 0;
            for j = 1:N(dir3)
                for i = 1:N(dir2)
                    m = m + 1;
                    pts(dir1,m) = end_pt(1);
                    pts(dir2,m) = x1D(i,N(dir2));
                    pts(dir3,m) = x1D(j,N(dir3));
                    w2 = w1D(j,N(dir2)) * ( N(dir2) > 1 ) + ( N(dir2) <= 1 );
                    w3 = w1D(j,N(dir3)) * ( N(dir3) > 1 ) + ( N(dir3) <= 1 );
                    wts(m) = w2 * w3;
                end
            end
        end
        function [pts,wts] = volm_quad_3D(N,x1D,w1D)
            pts = zeros( 3, prod(N) );
            wts = zeros( prod(N), 1 );
            m = 0;
            for k = 1:N(3)
                for j = 1:N(2)
                    for i = 1:N(1)
                        m = m+1;
                        pts(1,m) = x1D(i,N(1));
                        pts(2,m) = x1D(j,N(2));
                        pts(3,m) = x1D(k,N(3));
                        w1 = w1D(i,N(1)) * ( N(1) > 1 ) + ( N(1) <= 1 );
                        w2 = w1D(j,N(2)) * ( N(2) > 1 ) + ( N(2) <= 1 );
                        w3 = w1D(k,N(3)) * ( N(3) > 1 ) + ( N(3) <= 1 );
                        wts(m) = w1 * w2 * w3;
                    end
                end
            end
        end
    end

end