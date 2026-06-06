function [delta,branch] = vinokur_get_delta_two_sided(s0,s1,refine)
B   = sqrt(s0*s1);
tol = 0.001;
if abs(B-1) < tol
    delta  = 0;
    branch = 0;
elseif B > 1
    delta = vinokur_sinh_function( B );
    if (refine)
        delta = fzero(@(z)sinh(z)./z - B,delta);
    end
    branch = 1;
elseif B < 1
    delta = vinokur_sine_function( B );
    if (refine)
        delta = fzero(@(z)sin(z)./z - B,delta);
    end
    branch = 2;
end