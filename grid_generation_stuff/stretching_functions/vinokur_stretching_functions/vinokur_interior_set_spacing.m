function si = vinokur_interior_set_spacing(N,ti,di,refine)
dt = 1/(N-1);
s = di/dt;
% s = 1/s;
si = fzero(@(s)obj_fun(N,ti,s,di,refine),s);
end

function f = obj_fun(N,ti,si,di,refine)
    dt = 1/(N-1);
    [d2,~,~] = vinokur_interior(ti+dt,ti,si,refine);
    [d1,~,~] = vinokur_interior(ti,ti,si,refine);
    f = di - (d2-d1);
end