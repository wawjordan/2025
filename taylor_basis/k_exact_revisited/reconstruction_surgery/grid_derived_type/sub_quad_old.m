classdef sub_quad_old
properties
    pts (:,:) {mustBeNumeric}
    wts (:,1) {mustBeNumeric}
end
methods
    function this = sub_quad_old()
    end
    function this = allocate_sub_quad(this,n_dim, n_pts)
        this.pts = zeros(n_dim,n_pts);
        this.wts = zeros(n_pts,1);
    end
end
end