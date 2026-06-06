function [s0,s1] = vinokur_two_sided_set_one_spacing(N,d0,s1,refine)
dt = 1/(N-1);
s =d0/dt;
s0 = fzero(@(s)obj_fun(N,s,s1,d0,refine),s);
end

function f = obj_fun(N,s0,s1,d0,refine)
    dt = 1/(N-1);
    [d,~,~] = vinokur_two_sided(dt,s0,s1,refine);
    f = d-d0;
end