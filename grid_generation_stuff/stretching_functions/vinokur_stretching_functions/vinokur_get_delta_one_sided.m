function [delta,branch] = vinokur_get_delta_one_sided(s0,refine)
tol = 0.001;
if abs(s0-1) < tol
    delta  = 0;
    branch = 0;
elseif s0 > 1
    delta = 0.5*vinokur_sinh_function( s0 );
    if (refine)
        delta = 0.5*fzero(@(z)sinh(z)./z - s0,2*delta);
    end
    branch = 1;
elseif s0 < 1
    delta = 0.5*vinokur_sine_function( s0 );
    if (refine)
        delta = 0.5*fzero(@(z)sin(z)./z - s0,2*delta);
    end
    branch = 2;
end
end