function s0 = vinokur_one_sided_set_spacing(N,d0,refine)
dt = 1/(N-1);
s =d0/dt;
s0 = fzero(@(s)obj_fun(N,s,d0,refine),s);
end

function f = obj_fun(N,s0,d0,refine)
    [d,~,~] = vinokur_one_sided(1/(N-1),s0,refine);
    f = d-d0;
end