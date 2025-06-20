classdef face_vec_t_old
    properties
        n  (1,1) {mustBeInteger}
        v  (:,:) {mustBeNumeric}
        
    end

    methods
        function this = face_vec_t_old()
        end
        function this = set_face_normals_2D(this,X1,X2,X3,xq,dir,interpolant)
            sz = size(xq,2);
            this.n = sz;
            this.v = zeros(3,sz);
            for i = 1:sz
                Ja = interpolant.metrics(X1,X2,X3,xq(:,i));
                Ja(3,:) = 0;
                this.v(:,i) = Ja(:,dir)/norm( Ja(:,dir), 2 );
                kappa = interpolant.normal_curvature(X1,X2,X3,xq(:,i),dir,[0,1,0]);
            end
            
        end
        function this = set_face_normals_3D(this,X1,X2,X3,xq,dir,interpolant)
            sz = size(xq,2);
            this.n = sz;
            this.v = zeros(3,sz);
            for i = 1:sz
                Ja = interpolant.metrics(X1,X2,X3,xq(:,i));
                this.v(:,i) = Ja(:,dir)/norm( Ja(:,dir), 2 );
            end
        end
    end % methods
end