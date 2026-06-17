function val = c_grid_radial_spacing(xi,eta,wall_spacing,wake_spacing,boundary_distance,imax,jmax,num_wake_pts)
    
    b = (num_wake_pts-1)/(imax-1);
    a = b;% - 3/(imax-1);
    c = (imax-num_wake_pts)/(imax-1);
    d = c;% + 3/(imax-1);
    dh = wall_spacing;
    h  = wake_spacing;
    % [x,r,dx1,fh] = geomspace_4( N, xmin, xmax, dx0 )
    eta_min = @(xi) smooth_expand_outside(xi,jmax,jmax,a,b,c,d,dh,dh,dh,h,h,true);
    eta_min1 = @(xi) smooth_expand_outside2(xi,jmax,jmax,a,b,c,d,dh,dh,dh,h,h);
    % eta_max = @(xi) boundary_distance*ones(size(xi));
    r = @(xi) arrayfun(@(xi)rval(xi,jmax,boundary_distance,eta_min),xi);
    val = 0 + eta_min(xi).*( (r(xi).^(jmax-1)).^(eta) - 1 )./(r(xi)-1);
    % xi_min = @(eta) 0 + h*( (r^(jmax-1)).^(eta) - 1 )./(r-1);
    % xi_max = @(eta) 0 + h*( (r^(jmax-1)).^(eta) - 1 )./(r-1);

    % val = transfinite_surface_interpolate( xi, eta, ...
    %                                       @(eta) xi_min(eta), ...
    %                                       @(eta) xi_max(eta), ...
    %                                       @(xi) eta_min(xi),  ...
    %                                       @(xi) eta_max(xi) );
end

function r = rval(xi,jmax,boundary_distance,hfun)
dx0 = hfun(xi);
[~,r,~] = my_geomspace(jmax,hfun(0),xmax=boundary_distance,dx0=dx0);
end